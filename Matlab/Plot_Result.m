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
experflag=1;



SimMusclename=["bflh_r","bfsh_r","gaslat_r","gasmed_r","grac_r","recfem_r","sart_r","semimem_r","semiten_r","soleus_r","tfl_r","vasint_r","vaslat_r","vasmed_r"];
ExpMuscle=["RBICF","RSEMT","RMGAS","RRECF","RVASL","RVASM"];


%% Finding maximum EMG
   
load ([folder '\Data\FinalData_Seperated.mat']);
%%
TrialCounter=0;
if experflag
for A=1:length(AnalyzeMethod)
    for m=1:length(Modelname)
        results_folder = append(folder,"Result\",Modelname(m),"\");
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        for T1=1:length(Terials1)
            for T2=1:length(Terials2)
                TrialCounter=TrialCounter+1;
                filename=append(Terials1(T1),"_",Terials2(T2));

                if muscleplot
                figure
                x=RefT;
                t = tiledlayout(4,2);
                nexttile;
                title(t,[filename])
                ExpMotion=ResultData.(filename).NormalMotion; %scaling
                ExpForce=ResultData.(filename).NormalExForce;
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
                    
                    ExpMucs=ResultData.(filename).NormalExpEMG.(ExpMuscle(Flexmus)); %scaling
                    SimMusc=ResultData.(filename).NormalSimEMG.(SimMusclename(Flexmus));
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
end
for T1=1:length(Terials1)
    for T2=1:length(Terials2)
        TrialCounter=TrialCounter+1;
        filename=append(Terials1(T1),"_",Terials2(T2));
        for Flmus=1:length(ExpMuscle)
            ResultData.(filename).('ExpEMG').(ExpMuscle(Flmus));
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
Y=ResultData.MaxEMG.data([kk:kk+3],:);
X=[10 30 60 90];
plot(X,Y,'-*');
legend(ExpMuscle)
title(a(qq));
end
e = tiledlayout(1,2);
%%
ratio=ResultData.MaxEMG.data([1:8],:)./ResultData.MaxEMG.data([9:16],:);
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
Y=ResultData.MaxActivation.data([kk:kk+3],[2,10,5,6,14,15]);
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
Y=ResultData.MaxActuateForce.data([kk:kk+3],[2,10,5,6,14,15]);
bar(X,Y','stacked')
legend(ExpMuscle)
title(a(qq));
end

