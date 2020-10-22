% clear all;
% close all;
clc;
readflage=0;
% folder =[fileparts(mfilename('fullpath')) '\Data\'];
folder= 'C:\MyCloud\OneDriveUcf\Real\Simulation\Data\S1\RawData\';
fname = 'P005_T001_RKnee_';
Terials1=["Fl","Ex"];
Terials2=["IsoK60","IsoK120","IsoK180","IsoK240","IsoM10","IsoM30","IsoM60","IsoM90"];
% Terials2=["IsoM30","IsoM60","IsoM90"];
Fdata=[];
Gdata=[0];
k=0;
%
if readflage 
    for T1=1:length(Terials1)
        k=k+1;
        for T2=1:length(Terials2)
            
            Namedr(k)=append(Terials1(T1),"_",Terials2(T2));
            Datadr=append(folder,fname,Terials1(T1),"_",Terials2(T2),".csv");
            data=importdata(Datadr);
            [rf,cf]=find(strncmp(data.textdata,'Biodex',6));
            [rg,cg]=find(strncmp(data.textdata,'Gn',2));
            [rt,cT]=find(strncmp(data.textdata,'Biodex',6)|strncmp(data.textdata,'Gn',2));
            sr = data.data(3,cf-1)-data.data(2,cf-1);%finds sampling rates force
            srg =data.data(3,cT-1)-data.data(2,cT-1);%finds sampling rates Gonio
            for ii=1:length(cT)
                kk=1;
                for jj=1:length(data.data)%# of data points in a given trial
                    if ~isnan(data.data(jj,cT(ii)))&&data.data(jj,cT(ii))~=0
                        kk=jj;
                    end
                end
                t=srg(ii)*(kk-1); %time duration of trial based on # of data points
                y=interp1(0:srg(ii):t,data.data(1:kk,cT(ii)),linspace(0,t,t/(sr(1))+1)); %interpolates data to match sampling rate to force data
                b=y';
                ends(ii)=length(b);
                if (size(Gdata(:,1)) == 1) %recombines data into a matrix padded with NaN
                    Gdata = [[0:sr(1):t]' b];
                elseif (length(b) <  length(Gdata(:,1)))
                    Gdata = [Gdata [b; NaN(length(Gdata(:,1))-length(b),1)]];
                elseif (length(b) == length(Gdata(:,1)))
                    Gdata = [Gdata b];
                elseif (length(b) >  length(Gdata(:,1)))
                    Gdata = [[Gdata; NaN(length(b)-length(Gdata(:,1)),length(Gdata(1,:)))] b];
                end
                
            end
            FinalData.(Namedr(k)).data=Gdata;
            FinalData.(Namedr(k)).colheaders=["time" data.textdata(cT)];
            clear Gdata
            Gdata=[0];
            
        end
    end

save ([folder 'FinalData.mat'],'FinalData');

end
%%
load ([folder 'FinalData.mat']);
for T1=1:length(Terials1)
        k=k+1;
    for T2=1:length(Terials2)        
        Namedr(k)=append(Terials1(T1),"_",Terials2(T2));
        Data=FinalData.(Namedr(k)).data;
        HData=FinalData.(Namedr(k)).colheaders;
        [rg,cg]=find(strncmp(HData,'Gn',2));
        [rk,ck]=find(strncmp(HData,'Gn K',4));
        [rh,ch]=find(strncmp(HData,'Gn H',4));
        [rb,cb]=find(strncmp(HData,'Biodex',6));
        [r,c]=size(Data);
%% Process on Motion Data
        Gon=Data(:,cg);
        CalGon=-0.0058.*Gon.^2-1.699.*Gon;
        Data(:,cg)=CalGon.*pi()./180;
%% Save Motion
        fnames=['P005_T001_Motion_',char(Namedr{k}),'.mot'];
        fid=fopen([folder char(fnames)], "w");
        if fid < 0
            fprintf('\nERROR: %s could not be opened for writing...\n\n', name);
            return
        end
        fprintf(fid,[char(fnames) '\nversion=1\nnRows=%d\nnColumns=%d\nInDegrees=no\nendheader\n'],r,c);
        fprintf(fid,'time\tpelvis_tilt\tpelvis_tx\tpelvis_ty\thip_flexion_r\tknee_angle_r\tankle_angle_r\n');
        for i = 1:r
            fprintf(fid,'%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t',Data(i,1),0,0.055,1.059,Data(i,ch(2)),Data(i,ck(2)),0);
            fprintf(fid, '\n');
        end
            fclose(fid);
            fprintf('Saved %s\n', fnames)
%% Process Force
            A=[];
            x=1*Data(:,cb(1)); %data of a trial
            if x<0.176
                Mb=(-142.25).*x+24.8;
            else
                Mb=142.07.*x-25.32;
            end
%% Save Force
            fnames=['P005_T001_Torque_',char(Namedr{k}),'.mot'];
            fid=fopen([folder char(fnames)], "w");
            if fid < 0
                fprintf('\nERROR: %s could not be opened for writing...\n\n', name);
            return
            end
            fprintf(fid,[char(fnames) '\nversion=1\nnRows=%d\nnColumns=%d\nInDegrees=no\nendheader\n'],r,c);
            fprintf(fid,['time\treaction_force_vx\treaction_force_vy\treaction_force_vz\treaction_force_px\treaction_force_py\treaction_force_pz\treaction_torque_x\treaction_torque_y\treaction_torque_z\n']);
                for i = 1:r
                    fprintf(fid,'%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t',Data(i,1),0,0,0,0,0,0,0,0,Mb(i));
                    fprintf(fid, '\n');
                end
            fclose(fid);
            fprintf('Saved %s\n', fnames)
    end
end
