clear all
folder = 'C:\MyCloud\OneDriveUcf\Real\Simulation\Source\T001\';
AnalyzeMethod=["SOP","CMC"];
Modelname={'Rajagopal'};
Terials1=["Fl","Ex"];
Terials2=["IsoK60","IsoK120","IsoK180","IsoK240","IsoM10","IsoM30","IsoM60","IsoM90"];
Terials3=["iter1","iter2","iter3"];
TimeT=zeros(4,2);
delimiterIn = ',';
M_ThresholdMin=10*3.14/180;
M_ThresholdMax=90*3.14/180;
ForceRatio=0.4;
TrialCounter=0;
readflage=1;
muscleplot=0;
experflag=0;



SimMusclename=["bflh_r","bfsh_r","gaslat_r","gasmed_r","grac_r","recfem_r","sart_r","semimem_r","semiten_r","soleus_r","tfl_r","vasint_r","vaslat_r","vasmed_r"];
ExpMuscle=["RBICF","RSEMT","RMGAS","RRECF","RVASL","RVASM"];

if readflage
    for A=2:length(AnalyzeMethod)
        for m=1:length(Modelname)
            results_folder = append(folder,"Result\",Modelname(m),"\");
            for T1=1:length(Terials1)
                for T2=1:length(Terials2)
                    TrialCounter=TrialCounter+1;
                    filename=append(Terials1(T1),"_",Terials2(T2));
                    EMGFile=append(folder,"Data\P005_T001_Rknee_",filename,"_EMG.csv");
                    ExForcefile=append(folder,"Data\P005_T001_Rknee_",filename,"_Torque.mot");
                    IkFile=append(folder,"Data\P005_T001_Rknee_",filename,"_Motion.mot");
                    FTable=importdata(ExForcefile);
                    MTable=importdata(IkFile);
                    EMGtable=importdata(EMGFile,delimiterIn);
                    ExpRtime=FTable.data(:,1);
                    T01ResultData.(filename).('ExForce').('full')=FTable.data(:,[1,10]);
                    T01ResultData.(filename).('Motion').('full')=MTable.data(:,[1,6]);
                    T01ResultData.(filename).('ExEMG').('full')=EMGtable.data;
                    %% Finding time of iterative
                    Event=EventDetection(filename,FTable,ForceRatio,MTable,[M_ThresholdMin M_ThresholdMax]);
                    Stime=Event(:,1);
                    Etime=Event(:,2);
                    if length(Stime)~=3||length(Etime)~=3
                        fprintf('\nERROR: %s Wrong trail ...\n\n', filename);
                    end
                    MeanEmg=[];
                    for itr=1:length(Stime)
                        results_folder2=append(results_folder,AnalyzeMethod(A),"\",filename,"\",Terials3(itr),"\");
                        Expindx=find(ExpRtime>=Stime(itr)&ExpRtime<=Etime(itr));
                        T01ResultData.(filename).('time').(Terials3(itr))=MTable.data(Expindx,1);
                        T01ResultData.(filename).('Motion').(Terials3(itr))=MTable.data(Expindx,6);
                        T01ResultData.(filename).('ExForce').(Terials3(itr))=FTable.data(Expindx,10);
                        ActuateForce=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_Actuation_force.sto'));
                        MuscleActivation=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_controls.sto'));
                        Simtime=ActuateForce.data(1:(end-1),1); %Simulation time
                        Exptime=T01ResultData.(filename).('time').(Terials3(itr));
                        RefT=0:.1:100;
                        SimNormalT=((Simtime-Simtime(1))./(Simtime(end)-Simtime(1))).*100;
                        ExpNormalT=((Exptime-Exptime(1))./(Exptime(end)-Exptime(1))).*100;
                        y=interp1(SimNormalT,MuscleActivation.data(1:(end-1),:),RefT','linear','extrap');
                        y1=interp1(SimNormalT,ActuateForce.data(1:(end-1),:),RefT','linear','extrap');
                        y4=interp1(ExpNormalT,T01ResultData.(filename).ExForce.(Terials3(itr)),RefT','linear','extrap');
                        T01ResultData.(filename).('NormalExForce')(:,itr)=y4;
                        y4=interp1(ExpNormalT,T01ResultData.(filename).Motion.(Terials3(itr)),RefT','linear','extrap');
                        T01ResultData.(filename).('NormalMotion')(:,itr)=y4;
                        for Flexmus=1:length(ExpMuscle)
                            T01ResultData.(filename).('ExpEMG').(ExpMuscle(Flexmus)).(Terials3(itr))=EMGtable.data(Expindx,strncmp(EMGtable.colheaders,ExpMuscle(Flexmus),5));
                            MeanEmg(itr,Flexmus)=mean(EMGtable.data(Expindx,strncmp(EMGtable.colheaders,ExpMuscle(Flexmus),5)));
                            y3=interp1(ExpNormalT,T01ResultData.(filename).ExpEMG.(ExpMuscle(Flexmus)).(Terials3(itr)),RefT','linear','extrap');
                            T01ResultData.(filename).('NormalExpEMG').(ExpMuscle(Flexmus))(:,itr)=y3;
                        end
                        for Flmus=1:length(SimMusclename)
                            T01ResultData.(filename).('NormalSimEMG').(SimMusclename(Flmus))(:,itr)=y(:,strncmp(MuscleActivation.colheaders,SimMusclename(Flmus),5));
                            T01ResultData.(filename).('NormalSimForce').(SimMusclename(Flmus))(:,itr)=y1(:,strncmp(ActuateForce.colheaders,SimMusclename(Flmus),5));
                        end
                        
                        
                    end
                    MaxEMG(TrialCounter,:)=mean(MeanEmg);
                    
                end
            end
        end
    end
    %% Finding maximum EMG
    T01ResultData.('MaxEMG').colheaders= ExpMuscle;
    T01ResultData.('MaxEMG').data=MaxEMG;
    save ([folder '\Data\T01ResultData.mat'],'T01ResultData');
end
load ([folder '\Data\FinalData_Seperated.mat']);
%%
TrialCounter=0;
if experflag
for A=2:length(AnalyzeMethod)
    for m=1:length(Modelname)
        results_folder = append(folder,"Result\",Modelname(m),"\");
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        for T1=1:length(Terials1)
            for T2=1:length(Terials2)
                TrialCounter=TrialCounter+1;
                filename=append(Terials1(T1),"_",Terials2(T2));
%                 for itr=1:length(fieldnames(ResultData.(filename).time))
%                     results_folder2=append(results_folder,AnalyzeMethod(A),"\",filename,"\",Terials3(itr),"\");
%                     ActuateForce=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_Actuation_force.sto'));
%                     MuscleActivation=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_controls.sto'));
%                     Simtime=ActuateForce.data(1:(end-1),1); %Simulation time
%                     ExpRtime=ResultData.(filename).('time').(Terials3(itr));
%                     
%                     MaxActivation(TrialCounter,:)=mean(MuscleActivation.data);
%                     MaxActuateForce(TrialCounter,:)=mean(ActuateForce.data);
%                     
%                     
%                      RefT=0:.1:100;
%                      SimNormalT=((Simtime-Simtime(1))./(Simtime(end)-Simtime(1))).*100;
%                      
%                      y=interp1(SimNormalT,MuscleActivation.data(1:(end-1),:),RefT','linear','extrap');
%                      y1=interp1(SimNormalT,ActuateForce.data(1:(end-1),:),RefT','linear','extrap');
%                      for Flmus=1:length(SimMusclename)
%                          ResultData.(filename).('NormalSimEMG').(SimMusclename(Flmus))(:,itr)=y(:,strncmp(MuscleActivation.colheaders,SimMusclename(Flmus),5));
%                          ResultData.(filename).('NormalSimForce').(SimMusclename(Flmus))(:,itr)=y1(:,strncmp(ActuateForce.colheaders,SimMusclename(Flmus),5));
%                      end
%                      ExpNormalT=((ExpRtime-ExpRtime(1))./(ExpRtime(end)-ExpRtime(1))).*100;
%                    
%                      for Flexmus=1:length(ExpMuscle)
%                         y3=interp1(ExpNormalT,ResultData.(filename).ExpEMG.(ExpMuscle(Flexmus)).(Terials3(itr)),RefT','linear','extrap');
%                         ResultData.(filename).('NormalExpEMG').(ExpMuscle(Flexmus))(:,itr)=y3;
%                      end
%                        y4=interp1(ExpNormalT,ResultData.(filename).ExForce.(Terials3(itr)),RefT','linear','extrap');
%                        ResultData.(filename).('NormalExForce')(:,itr)=y4;
%                        y4=interp1(ExpNormalT,ResultData.(filename).Motion.(Terials3(itr)),RefT','linear','extrap');
%                        ResultData.(filename).('NormalMotion')(:,itr)=y4;
%                 end
                if muscleplot
                figure
                x=RefT;
                t = tiledlayout(4,2);
                nexttile;
                title(t,[filename])
                ExpMotion=T01ResultData.(filename).NormalMotion; %scaling
                ExpForce=T01ResultData.(filename).NormalExForce;
                RA1=mean(ExpMotion,2)-std(ExpMotion,0,2);
                RA2=mean(ExpMotion,2)+std(ExpMotion,0,2);
                RA3=mean(ExpMotion,2);
                                
                LA1=mean(ExpForce,2)-std(ExpForce,0,2);
                LA2=mean(ExpForce,2)+std(ExpForce,0,2);
                LA3=mean(ExpForce,2);
                x2 = [x, fliplr(x)];
                
                    RABetween = [RA1', fliplr(RA2')];
                    colororder([0.4940 0.1840 0.5560;0.4660 0.6740 0.1880])
                    yyaxis left
                    fill(x2, RABetween,[0.6 1 0.6]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                    hold on
                    plot(x,RA3,'color',[[0.4660 0.6740 0.1880]],'LineStyle','-');
                    ylabel('Force(N)');
%                     hold off
                    
                    %                 hold on
                    yyaxis right
                    LABetween = [LA1', fliplr(LA2')];
                    fill(x2, LABetween,[1 0.5 1]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                    plot(x,LA3,'color',[0.4940 0.1840 0.5560],'LineStyle','-');
                    xlabel('Normalized Time (%)');
                    ylabel('Angle(degree)');
                    hold off
                    
                for Flexmus=1:length(ExpMuscle)
                    
                    nexttile
                    
                    ExpMucs=T01ResultData.(filename).NormalExpEMG.(ExpMuscle(Flexmus)); %scaling
                    SimMusc=T01ResultData.(filename).NormalSimEMG.(SimMusclename(Flexmus));
                    %                 if strncmp(Terials1(T1),"Fl",2)
                    RA1=mean(ExpMucs,2)-std(ExpMucs,0,2);
                    RA2=mean(ExpMucs,2)+std(ExpMucs,0,2);
                    RA3=mean(ExpMucs,2);
                    
                    
                    LA1=mean(SimMusc,2)-std(SimMusc,0,2);
                    LA2=mean(SimMusc,2)+std(SimMusc,0,2);
                    LA3=mean(SimMusc,2);
                    
                    x2 = [x, fliplr(x)];
                    RABetween = [RA1', fliplr(RA2')];
                    yyaxis right
                    fill(x2, RABetween,[0.6 1 0.6]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                    hold on
                    plot(x,RA3,'color',[0.4660 0.6740 0.1880],'LineStyle','-');
%                     ylabel('Experiment activation');
                    
                    hold on
                    yyaxis left
                    LABetween = [LA1', fliplr(LA2')];
                    fill(x2, LABetween,[1 0.5 1]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                    plot(x,LA3,'color',[0.4940 0.1840 0.5560],'LineStyle','-');
                    title([SimMusclename(Flexmus)])
                    xlabel('Normalized Time (%)');
                    ylabel('activation ');
%                     legend('SD of Exp','SD of Sim','Mean of Exp','Mean of Sim');
                    hold off
                end
                end
%                 figure
%                 x=RefT;

                
            end
        end
    end
end
T01ResultData.('MaxActivation').colheaders= SimMusclename;
T01ResultData.('MaxActivation').data=MaxActivation;
T01ResultData.('MaxActuateForce').colheaders= SimMusclename;
T01ResultData.('MaxActuateForce').data=MaxActuateForce;
end
for T1=1:length(Terials1)
    for T2=1:length(Terials2)
        TrialCounter=TrialCounter+1;
        filename=append(Terials1(T1),"_",Terials2(T2));
        for Flmus=1:length(ExpMuscle)
            T01ResultData.(filename).('ExpEMG').(ExpMuscle(Flmus));
            MeanEmg(k,Flmus)=mean(EMGtable.data(Expindx,strncmp(EMGtable.colheaders,ExpMuscle(Flmus),5)));
        end
    end
end
 a=["Fl_IsoK","Fl_IsoM","Ex_IsoK","Ex_IsoM"'];
 e = tiledlayout(2,2);
 qq=0;
for kk=1:4:16
qq=qq+1;
nexttile(e);
Y=T01ResultData.MaxEMG.data([kk:kk+3],:);
X=[10 30 60 90];
plot(X,Y,'-*');
legend(ExpMuscle)
title(a(qq));
end
e = tiledlayout(1,2);
%%
ratio=T01ResultData.MaxEMG.data([1:8],:)./T01ResultData.MaxEMG.data([9:16],:);
qq=0;
for kk=1:4:8
qq=qq+1;
nexttile(e);
% Y=ResultData.EMG.data([kk:kk+3],:);
X=[10 30 60 90];
plot(X,ratio,'-*');
legend(ExpMuscle)
title(["IsoK","IsoM"]);
end
%%
for kk=1:4:16
qq=qq+1;
nexttile(e);
Y=T01ResultData.MaxActivation.data([kk:kk+3],[2,10,5,6,14,15]);
X=[10 30 60 90];
plot(X,Y,'-*');
legend(ExpMuscle)
title(a(qq));
end
X=[10 30 60 90];
 a=["Fl_IsoK","Fl_IsoM","Ex_IsoK","Ex_IsoM"'];
e = tiledlayout(2,2);
 qq=0;
for kk=1:4:16
qq=qq+1;
nexttile(e);
Y=T01ResultData.MaxActuateForce.data([kk:kk+3],[2,10,5,6,14,15]);
bar(X,Y','stacked')
legend(ExpMuscle)
title(a(qq));
end

