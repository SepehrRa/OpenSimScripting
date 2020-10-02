% clear all;
% close all;
clc;

folder =[fileparts(mfilename('fullpath')) '\Data\'];
fname = 'P005_T001_RKnee_';
Terials1=["EX";"FL"];
Terials2=["I60";"IK60"];
Fdata=[];
k=0;
%%
% for T1=1:length(Terials1)
%         k=k+1;
%     for T2=1:length(Terials2)
%         
%         Namedr(k)=append(Terials1(T1),"_",Terials2(T2));
%         Datadr=append(folder,fname,Terials1(T1),"_",Terials2(T2),".xlsx");
%         data=importdata(Datadr);
%         [rf,cf]=find(strncmp(data.textdata,'Biodex',6));
%         [rg,cg]=find(strncmp(data.textdata,'Gn',2));
%         [rt,cT]=find(strncmp(data.textdata,'Biodex',6)|strncmp(data.textdata,'Gn',2));
%         sr = data.data(3,cf-1)-data.data(2,cf-1);%finds sampling rates force 
%         srg =data.data(3,cT-1)-data.data(2,cT-1);%finds sampling rates Gonio
%             for ii=1:length(cT) 
%             kk=1;
%                 for jj=1:length(data.data)%# of data points in a given trial
%                     if ~isnan(data.data(jj,cT(ii)))&&data.data(jj,cT(ii))~=0
%                         kk=jj;
%                     end
%                 end 
%             t=srg(ii)*(kk-1); %time duration of trial based on # of data points
%             y=interp1(0:srg(ii):t,data.data(1:kk,cT(ii)),linspace(0,t,t/(sr(1))+1)); %interpolates data to match sampling rate to force data
%             b=y';
%             ends(ii)=length(b);
%                 if (size(Gdata(:,1)) == 1) %recombines data into a matrix padded with NaN
%                     Gdata = [[0:sr(1):t]' b];
%                 elseif (length(b) <  length(Gdata(:,1)))
%                     Gdata = [Gdata [b; NaN(length(Gdata(:,1))-length(b),1)]];
%                 elseif (length(b) == length(Gdata(:,1)))
%                     Gdata = [Gdata b];
%                 elseif (length(b) >  length(Gdata(:,1)))
%                     Gdata = [[Gdata; NaN(length(b)-length(Gdata(:,1)),length(Gdata(1,:)))] b];
%                 end
%                 
%             end
%             FinalData.(Namedr(k)).data=Gdata;
%             FinalData.(Namedr(k)).colheaders=["time" data.textdata(cT)];
%             clear Gdata
%             Gdata=[0];
%             
%             
% %             t=sr(1)*(kk-1);
% %             Gdata=[[0:sr(1):t]' Gdata];
% %             Fdata=[Fdata data.data(:,[cf(1)-1:cf(end)])];
%     end
% end
% save FinalData.mat FinalData;
% save Gdata.mat Gdata;
%%
load FinalData.mat
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
%       zer=mean(B(600:2600,1));
%       B=B-zer;
%       deg=mean(B(600:2600,4))./90;
%       B=B./deg;
        Data(:,cg)=Data(:,cg).*pi()./180;
%% Motion
        fnames=['P005_T001_Motion_',char(Namedr{k}),'.mot'];
        fid=fopen(char(fnames), "w");
        if fid < 0
            fprintf('\nERROR: %s could not be opened for writing...\n\n', name);
            return
        end
        fprintf(fid,[char(fnames) '\nversion=1\nnRows=%d\nnColumns=%d\nInDegrees=no\nendheader\n'],r,c);
        fprintf(fid,'time\tpelvis_tilt\tpelvis_tx\tpelvis_ty\thip_flexion_r\tknee_angle_r\tankle_angle_r\n');
        for i = 1:r
            fprintf(fid,'%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t',Data(i,1),0,0.055,1.059,Data(i,ch(1)),-Data(i,ck(1)),0);
            fprintf(fid, '\n');
        end
            fclose(fid);
            fprintf('Saved %s\n', fnames)
%% Force
            A=[];
%             l=.15;
%             mg=0;
            x=Data(:,cb(1)); %data of a trial
            Mb=6.0046.*x.^4+40.776.*x.^3+91.388.*x.^2-73.177.*x+9.5647; %finds moment about biodex arm
            %   Fr=Mb(1:rows(B))./l-mg.*sin((pi/2)-B(:,i))./2; %solves for reaction force
            %   A=[A Fr];
            fnames=['P005_T001_Torque_',char(Namedr{k}),'.mot'];
            fid=fopen(char(fnames), "w");
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
