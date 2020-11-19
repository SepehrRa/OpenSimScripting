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

SimMusclename=["bflh_r","semiten_r","gasmed_r","recfem_r","vaslat_r","vasmed_r","gaslat_r","gasmed_r","grac_r","sart_r","semimem_r","soleus_r","tfl_r","vasint_r","bfsh_r"];
ExpMuscle=["RBICF","RSEMT","RMGAS","RRECF","RVASL","RVASM"];
        %%
if readflage        
        for T1=2:length(Terials1)
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
               %% Finding time of iterative
               if strncmp(Terials2(T2),"IsoM",4)
                   [indx,c]=find(abs(FTable.data(:,10))>ForceRatio.*max(abs(FTable.data(:,10))));
                   Stime=FTable.data(indx([1;find(diff(indx)>10)+1]),1);
                   Etime=FTable.data(indx([(find(diff(indx)>10));end]),1);
         
                   %finding Maximum activation
                  
                   
               else
                   if strncmp(Terials1(T1),"Fl",2)
                       [indx,c]=find(MTable.data(:,6)>M_ThresholdMin & MTable.data(:,6)<M_ThresholdMax);
                       % & [0;diff(MTable.data(:,6))>0]
                       Sindx=indx([1;find(diff(indx)>10)+1]);
                       Eindx=indx([(find(diff(indx)>10));end]);
                       count=0;
                       Stime=[];
                       Etime=[];
                        for ww=1:length(Sindx)
                            if (MTable.data(Sindx(ww),6)-MTable.data(Sindx(ww)-1,6))>0
                                count=count+1;
                                Stime=[Stime;FTable.data(Sindx(ww),1)];
                                Etime=[Etime;FTable.data(Eindx(ww),1)];
                            end
                        end                    
                    else
                        [indx,c]=find(MTable.data(:,6)>M_ThresholdMin & MTable.data(:,6)<M_ThresholdMax);
                        Sindx=indx([1;find(diff(indx)>10)+1]);
                        Eindx=indx([(find(diff(indx)>10));end]);
                        count=0;
                        Stime=[];
                        Etime=[];
                        for ww=1:length(Sindx)
                            if (MTable.data(Sindx(ww)+1,6)-MTable.data(Sindx(ww),6))<0
                                count=count+1;
                                Stime=[Stime;FTable.data(Sindx(ww),1)];
                                Etime=[Etime;FTable.data(Eindx(ww),1)];
                            end
                        end
                    end
               end
                
    
                if length(Stime)~=3||length(Etime)~=3
                    fprintf('\nERROR: %s Wrong trail ...\n\n', filename);
                end
                MeanEmg=[];
                for k=1:length(Stime)
                    Expindx=find(ExpRtime>=Stime(k)&ExpRtime<=Etime(k));

                    ResultData.(filename).('time').(Terials3(k))=MTable.data(Expindx,1);
                    ResultData.(filename).('Motion').(Terials3(k))=MTable.data(Expindx,6);
                    ResultData.(filename).('ExForce').(Terials3(k))=FTable.data(Expindx,10);
                   
                    for Flmus=1:length(ExpMuscle)
                        ResultData.(filename).('ExpEMG').(ExpMuscle(Flmus)).(Terials3(k))=EMGtable.data(Expindx,strncmp(EMGtable.colheaders,ExpMuscle(Flmus),5));
                        MeanEmg(k,Flmus)=mean(EMGtable.data(Expindx,strncmp(EMGtable.colheaders,ExpMuscle(Flmus),5)));
                    end
                    
                end
                MaxEMG(TrialCounter,:)=mean(MeanEmg);

            end
        end
         %% Finding maximum EMG
          ResultData.('MaxEMG').colheaders= ExpMuscle;
          ResultData.('MaxEMG').data=max(MaxEMG);
          save ([folder '\Data\FinalData_Seperated.mat'],'ResultData');
end          
load ([folder '\Data\FinalData_Seperated.mat']);
            %%
for A=2:length(AnalyzeMethod)
    for m=1:length(Modelname)
        results_folder = append(folder,"Result\",Modelname(m),"\");
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        for T1=1:length(Terials1)
            for T2=6:length(Terials2)
                filename=append(Terials1(T1),"_",Terials2(T2));
                for itr=1:length(fieldnames(ResultData.(filename).time))
                    results_folder2=append(results_folder,AnalyzeMethod(A),"\",filename,"\",Terials3(itr),"\");
                    ActuateForce=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_Actuation_force.sto'));
                    MuscleActivation=importdata(append(results_folder2,Modelname(m),'_',filename,'_',Terials3(itr),'_CMC','_controls.sto'));
                    Simtime=ActuateForce.data(1:(end-1),1); %Simulation time
%                     ExpRtime=FTable.data(:,1);      %Experimental time
                    ExpRtime=ResultData.(filename).('time').(Terials3(itr));
                    %                 
                    % %                 ResultData=cell(5,11,3);
                    %                 %MuscleForce=
                    %
%                     y=interp1(data.data(ww+1:kk,ii-1),data.data(ww+1:kk,ii),[ts:DStime:te],'linear','extrap');
%                     
%                     Simindx=find(Simtime>=Stime(k)&Simtime<=Etime(k));
%                     Expindx=find(ExpRtime>=Stime(k)&ExpRtime<=Etime(k));
%                   figure  
%                     plot(Simtime(Simindx),MuscleActivation.data(Simindx,2),ExpRtime(Expindx),normalize(EMGtable.data(Expindx,2),'range'))
%                     
%                     RefT=ExpRtime(Expindx);
                     RefT=0:.1:100;
                     SimNormalT=((Simtime-Simtime(1))./(Simtime(end)-Simtime(1))).*100;
                     
                     y=interp1(SimNormalT,MuscleActivation.data(1:(end-1),:),RefT','linear','extrap');
                     y1=interp1(SimNormalT,ActuateForce.data(1:(end-1),:),RefT','linear','extrap');
                     for Flmus=1:length(SimMusclename)
                         ResultData.(filename).('NormalSimEMG').(SimMusclename(Flmus))(:,itr)=y(:,strncmp(MuscleActivation.colheaders,SimMusclename(Flmus),5));
                         ResultData.(filename).('NormalSimForce').(SimMusclename(Flmus))(:,itr)=y1(:,strncmp(ActuateForce.colheaders,SimMusclename(Flmus),5));
                     end
                     ExpNormalT=((ExpRtime-ExpRtime(1))./(ExpRtime(end)-ExpRtime(1))).*100;
                   
                     for Flexmus=1:length(ExpMuscle)
                        y3=interp1(ExpNormalT,ResultData.(filename).ExpEMG.(ExpMuscle(Flexmus)).(Terials3(itr)),RefT','linear','extrap');
                        ResultData.(filename).('NormalExpEMG').(ExpMuscle(Flexmus))(:,itr)=y3;
                     end
                      
                end
                for Flexmus=1:length(ExpMuscle)
                figure
                x=RefT;
                ExpMucs=ResultData.(filename).NormalExpEMG.(ExpMuscle(Flexmus))/ResultData.MaxEMG.data(Flexmus);
                SimMusc=ResultData.(filename).NormalSimEMG.(SimMusclename(Flexmus));
%                 if strncmp(Terials1(T1),"Fl",2)
                RA1=mean(ExpMucs,2)-std(ExpMucs,0,2);
                RA2=mean(ExpMucs,2)+std(ExpMucs,0,2);
                RA3=mean(ExpMucs,2);
                
                
                LA1=mean(SimMusc,2)-std(SimMusc,0,2);
                LA2=mean(SimMusc,2)+std(SimMusc,0,2);
                LA3=mean(SimMusc,2);
                yyaxis left
                x2 = [x, fliplr(x)];
                RABetween = [RA1', fliplr(RA2')];
                fill(x2, RABetween,[0.6 1 0.6]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                LABetween = [LA1', fliplr(LA2')];
                plot(x,RA3,'color',[[0.4660 0.6740 0.1880]]);
                
               
%                 hold on
                yyaxis right
                fill(x2, LABetween,[1 0.5 1]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                title([filename SimMusclename(Flexmus)]);
                plot(x,LA3,'color',[0.4940 0.1840 0.5560]);
                xlabel('Normalized Time (%)');
                ylabel('activation ');
                legend('SD of Exp','SD of Sim','Mean of Exp','Mean of Sim');
                hold off
                end
%                 end
            end
        end
    end
end
