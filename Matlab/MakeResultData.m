clear all
folder = 'C:\MyCloud\OneDriveUcf\Real\Simulation\Source\T001\';
AnalyzeMethod=["SOP","CMC"];
Modelname={'Rajagopal'};
Terials1=["Fl","Ex"];
Terials2=["IsoK60","IsoK120","IsoK180","IsoK240","IsoM10","IsoM30","IsoM60","IsoM90"];
Terials3=["iter1","iter2","iter3"];
M_ThresholdMin=10*3.14/180;
M_ThresholdMax=90*3.14/180;
ForceRatio=0.4;
SimMusclename=["bflh_r","bfsh_r","gaslat_r","gasmed_r","grac_r","recfem_r","sart_r","semimem_r","semiten_r","soleus_r","tfl_r","vasint_r","vaslat_r","vasmed_r"];
ExpMuscle=["RBICF","RSEMT","RMGAS","RRECF","RVASL","RVASM"];
    for A=2:length(AnalyzeMethod)
        for m=1:length(Modelname)
            results_folder = append(folder,"Result\",Modelname(m),"\");
            for T1=1:length(Terials1)
                for T2=1:length(Terials2)
                    filename=append(Terials1(T1),"_",Terials2(T2));
                    EMGFile=append(folder,"Data\P005_T001_Rknee_",filename,"_EMG.csv");
                    ExForcefile=append(folder,"Data\P005_T001_Rknee_",filename,"_Torque.mot");
                    IkFile=append(folder,"Data\P005_T001_Rknee_",filename,"_Motion.mot");
                    FTable=importdata(ExForcefile);
                    MTable=importdata(IkFile);
                    EMGtable=importdata(EMGFile);
                    ExpRtime=FTable.data(:,1);
                    ResultData.(filename).('ExForce').('full')=FTable.data(:,[1,10]);
                    ResultData.(filename).('Motion').('full')=MTable.data(:,[1,6]);
                    ResultData.(filename).('ExEMG').('full')=EMGtable.data;
                    %% Finding time of iterative
                    Event=EventDetection(filename,FTable,ForceRatio,MTable,[M_ThresholdMin M_ThresholdMax]);
                    Stime=Event(:,1);
                    Etime=Event(:,2);
                    if length(Stime)~=3||length(Etime)~=3
                        fprintf('\nERROR: %s Wrong trail ...\n\n', filename);
                    end
                    for itr=1:length(Stime)
                        results_folder2=append(results_folder,AnalyzeMethod(A),"\",filename,"\",Terials3(itr),"\");
                        Expindx=find(ExpRtime>=Stime(itr)&ExpRtime<=Etime(itr));
                        ResultData.(filename).('time').(Terials3(itr))=MTable.data(Expindx,1);
                        ResultData.(filename).('Motion').(Terials3(itr))=MTable.data(Expindx,6);
                        ResultData.(filename).('ExForce').(Terials3(itr))=FTable.data(Expindx,10);
                        ActuateForce=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_Actuation_force.sto'));
                        MuscleActivation=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_controls.sto'));
                        Simtime=ActuateForce.data(1:(end-1),1); %Simulation time
                        Exptime=ResultData.(filename).('time').(Terials3(itr));
                        RefT=0:.1:100;
                        SimNormalT=((Simtime-Simtime(1))./(Simtime(end)-Simtime(1))).*100;
                        ExpNormalT=((Exptime-Exptime(1))./(Exptime(end)-Exptime(1))).*100;
                        y=interp1(SimNormalT,MuscleActivation.data(1:(end-1),:),RefT','linear','extrap');
                        y1=interp1(SimNormalT,ActuateForce.data(1:(end-1),:),RefT','linear','extrap');
                        
                        for Flmus=1:length(SimMusclename)
                            ResultData.(filename).('NormalSimEMG').(SimMusclename(Flmus))(:,itr)=y(:,strncmp(MuscleActivation.colheaders,SimMusclename(Flmus),5));
                            ResultData.(filename).('NormalSimForce').(SimMusclename(Flmus))(:,itr)=y1(:,strncmp(ActuateForce.colheaders,SimMusclename(Flmus),5));
                        end
                        y2=interp1(ExpNormalT,ResultData.(filename).ExForce.(Terials3(itr)),RefT','linear','extrap');
                        ResultData.(filename).('NormalExForce')(:,itr)=y2;
                        y3=interp1(ExpNormalT,ResultData.(filename).Motion.(Terials3(itr)),RefT','linear','extrap');
                        ResultData.(filename).('NormalMotion')(:,itr)=y3;
                        for Flexmus=1:length(ExpMuscle)
                            ResultData.(filename).('ExpEMG').(ExpMuscle(Flexmus)).(Terials3(itr))=EMGtable.data(Expindx,strncmp(EMGtable.colheaders,ExpMuscle(Flexmus),5));
                            MeanEmg(itr,Flexmus)=mean(EMGtable.data(Expindx,strncmp(EMGtable.colheaders,ExpMuscle(Flexmus),5)));
                            y4=interp1(ExpNormalT,ResultData.(filename).ExpEMG.(ExpMuscle(Flexmus)).(Terials3(itr)),RefT','linear','extrap');
                            ResultData.(filename).('NormalExpEMG').(ExpMuscle(Flexmus))(:,itr)=y4;
                        end

                        
                        
                    end                   
                end
            end
        end
    end
    save ([folder '\Data\T01ResultData.mat'],'ResultData');

