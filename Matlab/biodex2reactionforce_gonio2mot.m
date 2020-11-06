% clear all;
% close all;
clc;
 % Some times there is no need to import raw data because all data will
 % save in FinalDatafor first time. readflage=1 means import files again.
readflage=0; 
folder=uigetdir(); % get Data directory 
fname = 'P005_T001_RKnee_';
Terials1=["Ex","Fl"];
Terials2=["IsoK60","IsoK120","IsoK180","IsoK240","IsoM10","IsoM30","IsoM60","IsoM90"];
% Terials1=["Fl"];
% Terials2=["IsoK60"];
Fdata=[];
Gdata=[0];
k=0;
DStime=0.01; % desierd sampling time 
%
if readflage 
    for T1=1:length(Terials1)
        k=k+1;
        for T2=1:length(Terials2)
            
            Namedr(k)=append(Terials1(T1),'_',Terials2(T2));
            Datadr=append(folder,"\",fname,Terials1(T1),"_",Terials2(T2),".csv");
            data=importdata(Datadr);
            [rf,cf]=find(strncmp(data.textdata,'Biodex',6));
            [rg,cg]=find(strncmp(data.textdata,'Gn',2));
            [rt,cT]=find(strncmp(data.textdata,'Biodex',6)|strncmp(data.textdata,'Gn',2));
            sr = data.data(3,cf-1)-data.data(2,cf-1);%finds sampling rates force
            srg =data.data(3,cT-1)-data.data(2,cT-1);%finds sampling rates Gonio
            for ii=1:length(cT)
                kk=1;
                tflage=0;
                ww=0;
                for jj=1:length(data.data)%# of data points in a given trial
                    if ~isnan(data.data(jj,cT(ii)))&&data.data(jj,cT(ii))~=0
                        kk=jj;  %finding zero data at the end of each chanel
                        tflage=1;                        
                    elseif tflage==0 
                        ww=jj;  %finding zero data at the begining 
                    end
                end

                if ii==1
                    ts=data.data(ww+1,1);
                    te=data.data(kk,1); % final time of first chanel to set as final time for every other channel.
                end 
                y=interp1(data.data(ww+1:kk,cT(ii)-1),data.data(ww+1:kk,cT(ii)),[ts:DStime:te]); %Interpolates data to match sampling time to desierd sampling time 
                b=y';
                ends(ii)=length(b);
                if (size(Gdata(:,1)) == 1) %recombines data into a matrix padded with NaN
                    Gdata = [[data.data(ww+1,1):DStime:data.data(kk,1)]' b];
                else 
                    Gdata = [Gdata b];
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
        %find Knee Goniometer
        [rk,ck]=find(strncmp(HData,'Gn K',4));
        %find Hip goniometer
        [rh,ch]=find(strncmp(HData,'Gn H',4));
        [rb,cb]=find(strncmp(HData,'Biodex',6));
        [r,c]=size(Data);
%% Process on Motion Data
        Gon=Data(:,cg);
        CalGon=-0.0058.*Gon.^2-1.62.*Gon+1.14;  % Goniometer calibration, This equation would change base on ne calibration curve
        CalGon(CalGon<0)=0;
        Data(:,cg)=CalGon.*pi()./180;
%% Save Motion
        fnames=[fname,char(Namedr{k}),'_Motion.mot'];
        fid=fopen([folder '\' char(fnames)], "w");
        if fid < 0
            fprintf('\nERROR: %s could not be opened for writing...\n\n', fname);
            return
        end
        fprintf(fid,[char(fnames) '\nversion=1\nnRows=%d\nnColumns=%d\nInDegrees=no\nendheader\n'],r,7);
        fprintf(fid,'time\tpelvis_tilt\tpelvis_tx\tpelvis_ty\thip_flexion_r\tknee_angle_r\tankle_angle_r\n');
        for i = 1:r
            fprintf(fid,'%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t',Data(i,1),0,0.055,1.059,Data(i,ch(2)),Data(i,ck(2)),0);
            fprintf(fid, '\n');
        end
            fclose(fid);
            fprintf('Saved %s\n', [folder char(fnames)])
%% Process Force
            A=[];
            x=1*Data(:,cb(1)); %data of a trial
            Mb=-1.*(141.81.*x-25.047);
%% Save Force
            F_fnames=[fname,char(Namedr{k}),'_Torque.mot'];
            fid=fopen([folder '\' char(F_fnames)], "w");
            if fid < 0
                fprintf('\nERROR: %s could not be opened for writing...\n\n', fname);
            return
            end
            fprintf(fid,[char(F_fnames) '\nversion=1\nnRows=%d\nnColumns=%d\nInDegrees=no\nendheader\n'],r,10);
            fprintf(fid,['time\treaction_force_vx\treaction_force_vy\treaction_force_vz\treaction_force_px\treaction_force_py\treaction_force_pz\treaction_torque_x\treaction_torque_y\treaction_torque_z\n']);
                for i = 1:r
                    fprintf(fid,'%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t',Data(i,1),0,0,0,0,0,0,0,0,Mb(i));
                    fprintf(fid, '\n');
                end
            fclose(fid);
            fprintf('Saved %s\n', [folder char(F_fnames)])
    end
end
