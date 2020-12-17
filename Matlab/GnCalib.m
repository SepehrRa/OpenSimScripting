function [Ph,Pk,Pa]= GnCalib (folder)
close all
clear all
% folder = 'C:\MyCloud\OneDriveUcf\Real\Simulation\Source\T002\Data\';
Fname="P005_T002_GnCalib_";
Names=["AnkleDors20","AnkleDors10","Ankle0","AnklePlant10","AnklePlant30","AnklePlant50"...
    ,"Hip0","Hip30","Hip60","Hip90"...
    ,"Knee0","Knee30","Knee45","Knee60","Knee90"];
AnkelGnCalibdata=[];
HipGnCalibdata=[];
KneeGnCalibdata=[];

for i=1:length(Names)
    GnCaldatadir=importdata(append(folder,Fname,Names(i),".csv"));
    [Gr,Gc]=find(strncmp(GnCaldatadir.textdata,'Gn',2));
%     [Tr,Tc]=find(strncmp(GnCaldatadir.textdata,'X',1));
    Tc=[1,3];
    if length (GnCaldatadir.textdata)<4
        Tc=[1 1];
    end
    ChAindx=find(GnCaldatadir.data(:,Tc(1))<=2.5 & GnCaldatadir.data(:,Tc(1))>=1.5);
    ChA=mean(GnCaldatadir.data(ChAindx,Gc(1)));
    ChBindx=find(GnCaldatadir.data(:,Tc(2))<=2.5 & GnCaldatadir.data(:,Tc(2))>=1.5);
    ChB=mean(GnCaldatadir.data(ChBindx,Gc(2)));
    
    if contains(Names(i),"Ankle")
        AnkelGnCalibdata=[AnkelGnCalibdata;ChA ChB];
    elseif contains(Names(i),"Hip")
        HipGnCalibdata=[HipGnCalibdata;ChA ChB];
    elseif contains(Names(i),"Knee")
        KneeGnCalibdata=[KneeGnCalibdata;ChA ChB];
    end
    
end
x=[-20,-10,0,10,30,50];
Pa = polyfit(x,AnkelGnCalibdata(:,2)',2);
y1 = polyval(Pa,x);
plot(x,AnkelGnCalibdata(:,1),'*',x,y1,x,AnkelGnCalibdata(:,2),'*');
title('Ankle')
figure
x=[0,30,60,90];
Ph = polyfit(x,HipGnCalibdata(:,2)',2);
y1 = polyval(Ph,x);
plot(x,HipGnCalibdata(:,1),'*',x,y1,x,HipGnCalibdata(:,2),'*');
title('Hip')
figure
x=[0,30,45,60,90];
% x1 = linspace(0,4*pi)
Pk = polyfit(x,KneeGnCalibdata(:,2)',2);
y1 = polyval(Pk,x);
plot(x,KneeGnCalibdata(:,1),'*',x,y1,x,KneeGnCalibdata(:,2),'*',x,-sqrt(KneeGnCalibdata(:,2).^2+KneeGnCalibdata(:,1).^2),'*');
title('Knee')
end