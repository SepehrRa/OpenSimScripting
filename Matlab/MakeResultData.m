clear all
SubjectNumber='T002';
Project='P005';
folder=append('C:\MyCloud\OneDriveUcf\Real\Simulation\Source\',SubjectNumber);
psname=append(Project,'_',SubjectNumber);
legname='RKnee';
Subjectname =append(psname,"_",legname,"_");
AnalyzeMethod=["CMC"];
Modelname=["Rajagopal","Arnold"];
% Modelname=["Arnold"];
Terials1=["Fl"];
% Terials2=["IsoM10","IsoM30","IsoM45","IsoM60","IsoM90","IsoK60","IsoK120","IsoK180","IsoK240"];
Terials2=["IsoM10","IsoM30","IsoM45","IsoM60","IsoM90","IsoK60","IsoK120"];
Terials3=["iter1","iter2","iter3"];
results_folder = append(folder,"\Result");
load (append(results_folder,"\",psname,"_ResultData.mat"));
Coordinatesname=["knee_angle","ankle_angle","subtalar_angle"];
SimMusclename=["bflh_r","bfsh_r","gaslat_r","gasmed_r","recfem_r","sart_r","semimem_r","semiten_r","soleus_r","tfl_r","vasint_r","vaslat_r","vasmed_r"];
ExpMuscle=["RBICF","RSEMT","RMGAS","RRECF","RVASL","RVASM"];
CMCErr_Thereshold=2;
counter=0;
for m=1:length(Modelname)
    for A=1:length(AnalyzeMethod)   
        for T1=1:length(Terials1)
            for T2=1:length(Terials2)
                filename=append(Terials1(T1),"_",Terials2(T2));
                counter=counter+1;
                %                     EMGFile=append(folder,"\Data\",Subjectname,filename,"_EMG.csv");
                %                     ExForcefile=append(folder,"\Data\",Subjectname,filename,"_Torque.mot");
                %                     IkFile=append(folder,"\Data\",Subjectname,filename,"_Motion.mot");
                % Importing files
                %                     FTable=importdata(ExForcefile);
                %                     MTable=importdata(IkFile);
                %                     EMGtable=importdata(EMGFile);
                ExpRtime= ResultData.(filename).ExpForce.full(:,1);
                %                     ResultData.(filename).('ExpForce').('full')=FTable.data(:,[1,10]);
                %                     ResultData.(filename).('Motion').('full')=MTable.data(:,[1,6]);
                %                     ResultData.(filename).('ExpEMG').('full')=EMGtable.data;
                %% Finding events
                %                     Event=EventDetection(filename,FTable.data,ResultData.info.ForceRatio,MTable.data,[ResultData.info.M_ThresholdMin ResultData.info.M_ThresholdMax]);
                %                     Stime=Event(:,1);
                %                     Etime=Event(:,2);
                %                     %% Trail check
                %                     if length(Stime)~=3||length(Etime)~=3
                %                         fprintf('\nERROR: %s Wrong trail ...\n\n', filename);
                %                     end
                %% Strat reading Simulation files
                
                for itr=1:length(fieldnames(ResultData.(filename).time.Exp))
                    counter=counter+1;
                    Stime= ResultData.(filename).time.Exp.(Terials3(itr))(1);
                    Etime= ResultData.(filename).time.Exp.(Terials3(itr))(end);
                    %%% Making Directory
                    results_folder2=append(results_folder,"\",Modelname(m),"\",AnalyzeMethod(A),"\",filename,"\",Terials3(itr),"\");
                    Expindx=find(ExpRtime>=Stime&ExpRtime<=Etime);
                    %%% Reading Muscle Force, Muscle activation
                    ActuateForce=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_Actuation_force.sto'));
                    MuscleActivation=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_controls.sto'));
                    %%% Saving Experiment time when events happened
                    %                       ResultData.(filename).('time').Exp.(Terials3(itr))=MTable.data(Expindx,1);
                    %%% Saving Simulation time when events happened
                    ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('time').Sim.(Terials3(itr))=ActuateForce.data(1:(end-1),1);
                    %%% Saving Motion and Torque data when events happened
                    %                         ResultData.(filename).('Motion').(Terials3(itr))=MTable.data(Expindx,6);
                    %                         ResultData.(filename).('ExpForce').(Terials3(itr))=FTable.data(Expindx,10);
                    %%% Saving Motion and Torque data when events happened
                    Simtime=ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).time.Sim.(Terials3(itr)); %Simulation time
                    Exptime=ResultData.(filename).('time').Exp.(Terials3(itr));
                    %%% Defining reference time
                    RefT=0:.1:100;
                    %%% making time of each event in same scale of
                    %%% 0-100
                    SimNormalT=((Simtime-Simtime(1))./(Simtime(end)-Simtime(1))).*100;
                    %                       ExpNormalT=((Exptime-Exptime(1))./(Exptime(end)-Exptime(1))).*100;
                    %%% Interpolation to make time normalized data
                    y=interp1(SimNormalT,MuscleActivation.data(1:(end-1),:),RefT','linear','extrap');
                    y1=interp1(SimNormalT,ActuateForce.data(1:(end-1),:),RefT','linear','extrap');
                    %%% Saving time normalized of event's time from
                    %%% simualtion
                    ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('time').('NormalT')(:,itr)=y(:,1);
                    %%% Saving time normalized activation and force of
                    %%% simulation to each muscle
                    for Flmus=1:length(SimMusclename)
                        ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('SimActivation').(SimMusclename(Flmus)).(Terials3(itr))=MuscleActivation.data(:,strncmp(MuscleActivation.colheaders,SimMusclename(Flmus),5));
                        ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('NormalSimActivation').(SimMusclename(Flmus))(:,itr)=y(:,strncmp(MuscleActivation.colheaders,SimMusclename(Flmus),5));
                        
                        ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('SimMuscleForce').(SimMusclename(Flmus)).(Terials3(itr))=ActuateForce.data(:,strncmp(ActuateForce.colheaders,SimMusclename(Flmus),5));
                        ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('NormalSimMuscleForce').(SimMusclename(Flmus))(:,itr)=y1(:,strncmp(ActuateForce.colheaders,SimMusclename(Flmus),5));
                    end
                    for Coor=1:length(Coordinatesname)
                        CMCErr=ActuateForce.data(:,strcmpi(ActuateForce.colheaders,append(Coordinatesname(Coor),"_r_reserve")));
                        if mean(abs(CMCErr))> CMCErr_Thereshold
                             fprintf('\nWarning: %s has high residual number ...\n\n', append(filename,"_",Terials3(itr),"_",Coordinatesname(Coor)));
                        end
                        ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('Residual')((itr),Coor)=mean(CMCErr);
                        CMCErrTotal(counter,Coor)=mean(CMCErr);
                    end
                    %%% Error checking
                    
                    %%% Making time normalized external force
                    y2=interp1(Exptime,ResultData.(filename).ExpForce.(Terials3(itr)),ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).time.NormalT(:,itr),'linear','extrap');
                    ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('NormalExpForce')(:,itr)=y2;
                    %%% Making time normalized external motion
                    y3=interp1(Exptime,ResultData.(filename).Motion.(Terials3(itr)),ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('time').('NormalT')(:,itr),'linear','extrap');
                    ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('NormalMotion')(:,itr)=y3;
                    %%% Making time normalized experimental event time
                    %                       ENt=interp1(ExpNormalT,Exptime,RefT','linear','extrap');
                    %                       ResultData.(filename).('time').('ExpNormalT')(:,itr)=ENt;
                    %%% Making time normalized EMG data
                    y4=interp1(Exptime,ResultData.(filename).ExpEMG.full(Expindx,:),ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('time').('NormalT')(:,itr),'linear','extrap');
                    for Flexmus=1:length(ExpMuscle)
                        %%% Saving EMG data when events happened
%                         ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('ExpEMG').(ExpMuscle(Flexmus)).(Terials3(itr))=EMGtable.data(Expindx,strncmp(EMGtable.colheaders,ExpMuscle(Flexmus),5));
                        %%% Saving time normalized EMG of each muscle
                        ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).('NormalExpEMG').(ExpMuscle(Flexmus))(:,itr)=y4(:,strncmp(ResultData.(filename).ExpEMG.colheaders,ExpMuscle(Flexmus),5));
                    end
                    
                    
                    
                end
            end
        end
    end
end
save (append(results_folder,"\",psname,"_SimExpResultData.mat"),'ResultData');

