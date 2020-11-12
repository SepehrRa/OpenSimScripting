folder = 'C:\MyCloud\OneDriveUcf\Real\Simulation\Source\T001\';
AnalyzeMethod=["SOP","CMC"];
Modelname={'Rajagopal'};
Terials1=["Fl","Ex"];
Terials2=["IsoM10","IsoM30","IsoM60","IsoM90","IsoK60","IsoK120","IsoK180","IsoK240"];
TimeT=zeros(4,2);
delimiterIn = '\t';
M_ThresholdMin=10*3.14/180;
M_ThresholdMax=80*3.14/180;
% for A=1:length(AnalyzeMethod)
    for m=1:length(Modelname)
        results_folder = append(folder,"Result\",Modelname(m),"\");
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        CMCSetup=append('CMC\',Modelname(m),'\CMC_Setup_I0.xml');
        IDSetup=append('CMC\',Modelname(m),'\ID_Setup.xml');
        for T1=5:length(Terials1)
            for T2=1:length(Terials2)
                filename=append(Terials1(T1),"_",Terials2(T2));
                results_folder2=append(results_folder,AnalyzeMethod,"\",filename,"\");
                ExForcefile=append(folder,"Data\P005_T001_Rknee_",filename,"_Torque.mot");
                IkFile=append(folder,"Data\P005_T001_Rknee_",filename,"_Motion.mot");
                FTable=importdata(ExForcefile,delimiterIn);
                MTable=importdata(IkFile,delimiterIn);
                if strncmp(Terials2(T1),"IsoM",4)
                    [indx,c]=find(abs(FTable.data(:,10))>0.6*max(abs(FTable.data(:,10))));
                    Stime=FTable.data(indx([1;find(diff(indx)>10)]),1);
                    Etime=FTable.data(indx([(find(diff(indx)>10));end]),1);
                else
                    if strncmp(Terials2(T1),"Fl",2)
                        [indx,c]=find(MTable.data(:,6)>M_ThresholdMin||MTable.data(:,6)<M_ThresholdMax);
                        Stime=FTable.data(indx([1;find(diff(indx)>10)]),1);
                        Etime=FTable.data(indx([(find(diff(indx)>10));end]),1);
                    else
                    end
                end
            end
        end
    end
    