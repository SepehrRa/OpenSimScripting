folder = 'C:\MyCloud\OneDriveUcf\Real\Simulation\Source\T001\';
AnalyzeMethod=["SOP","CMC"];
Modelname={'Rajagopal'};
Terials1=["Fl","Ex"];
Terials2=["IsoM10","IsoM30","IsoM60","IsoM90","IsoK60","IsoK120","IsoK180","IsoK240"];
TimeT=zeros(4,2);
delimiterIn = '\t';
M_ThresholdMin=10*3.14/180;
M_ThresholdMax=80*3.14/180;
ForceRatio=0.5;
Musclename=["bflh_r","bfsh_r","gaslat_r","gasmed_r","grac_r","recfem_r","sart_r","semimem_r","semiten_r","vasint_r","vaslat_r","vasmed_r"];
for A=2:length(AnalyzeMethod)
    for m=1:length(Modelname)
        results_folder = append(folder,"Result\",Modelname(m),"\");
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        CMCSetup=append('CMC\',Modelname(m),'\CMC_Setup_I0.xml');
        IDSetup=append('CMC\',Modelname(m),'\ID_Setup.xml');
        for T1=2:length(Terials1)
            for T2=1:length(Terials2)
                filename=append(Terials1(T1),"_",Terials2(T2));
                results_folder2=append(results_folder,AnalyzeMethod(A),"\",filename,"\");
                ExForcefile=append(folder,"Data\P005_T001_Rknee_",filename,"_Torque.mot");
                IkFile=append(folder,"Data\P005_T001_Rknee_",filename,"_Motion.mot");
                EMGFile=append(folder,"Data\P005_T001_Rknee_",filename,"_EMG.mot");
                FTable=importdata(ExForcefile,delimiterIn);
                MTable=importdata(IkFile,delimiterIn);
                EMGable=importdata(EMGFile,delimiterIn);
                if strncmp(Terials2(T2),"IsoM",4)
                    [indx,c]=find(abs(FTable.data(:,10))>ForceRatio*max(abs(FTable.data(:,10))));
                    Stime=FTable.data(indx([1;find(diff(indx)>10)+1]),1);
                    Etime=FTable.data(indx([(find(diff(indx)>10));end]),1);
                else
                    if strncmp(Terials1(T1),"Fl",2)
                        [indx,c]=find(MTable.data(:,6)>M_ThresholdMin & MTable.data(:,6)<M_ThresholdMax & [0;diff(MTable.data(:,6))>0]);
                        Stime=FTable.data(indx([1;find(diff(indx)>10)+1]),1);
                        Etime=FTable.data(indx([(find(diff(indx)>10));end]),1);
                    else
                        [indx,c]=find(MTable.data(:,6)>M_ThresholdMin & MTable.data(:,6)<M_ThresholdMax & [0;diff(MTable.data(:,6))<0]);
                        Stime=FTable.data(indx([1;find(diff(indx)>10)+1]),1);
                        Etime=FTable.data(indx([(find(diff(indx)>10));end]),1);
                    end
                end
                if length(Stime)~=3||length(Etime)~=3
                    fprintf('\nERROR: %s Wrong trail ...\n\n', filename);
                end
                %%read CMC files
                ActuateForce=importdata(append(results_folder2,Modelname(m),"_",filename,"_",AnalyzeMethod(A),"_Actuation_force.sto"));
                MuscleActivation=importdata(append(results_folder2,Modelname(m),"_",filename,"_",AnalyzeMethod(A),"_controls.sto"));
                ARtime=ActuateForce.data(:,1);
                FRtime=FTable.data(:,1);
                RefT=0:.05:100;
                ResultData=cell(5,3);
                %MuscleForce=
                for k=1:length(Stime)
                    Aindx=find(ARtime>=Stime(k)&ARtime<=Etime(k));
                    Findx=find(FTable.data(:,1)>=Stime(k)&FTable.data(:,1)<=Etime(k));
                    NormalT=((ARtime(Aindx)-ARtime(Aindx(1)))./(ARtime(Aindx(end))-ARtime(Aindx(1)))).*100;
                    FNormalT=((FRtime(Findx)-FRtime(Findx(1)))./(FRtime(Findx(end))-FRtime(Findx(1)))).*100;
                    ResultData(1,k)={interp1(NormalT,ActuateForce.data(Aindx,2:length(Musclename)),RefT')};
                    ResultData(2,k)={interp1(NormalT,MuscleActivation.data(MuscleActivation,2:length(Musclename)),RefT')};
                    ResultData(3,k)={interp1(FNormalT,MTable.data(Findx,6),RefT')};
                    ResultData(4,k)={interp1(FNormalT,FTable.data(Findx,10),RefT')};
                    ResultData(5,k)={interp1(FNormalT,EMGable.data(Findx,2:end),RefT')};
                end
                
            end
        end
    end
end
