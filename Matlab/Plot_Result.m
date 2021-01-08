clear all
close all
% folder = 'C:\MyCloud\OneDriveUcf\Real\Simulation\Source\T002\';
AnalyzeMethod=["CMC"];
Modelname=["Rajagopal"];
Terials1=["Ex","Fl"];
Terials2=["IsoM10","IsoM30","IsoM45","IsoM60","IsoM90","IsoK60","IsoK120","IsoK180","IsoK240"];
Terials3=["iter1","iter2","iter3"];
TimeT=zeros(4,2);
delimiterIn = ',';
M_ThresholdMin=10*3.14/180;
M_ThresholdMax=90*3.14/180;
ForceRatio=0.8;
TrialCounter=0;
readflage=1;
muscleplot=1;
experflag=1;
Lcolorcode=[0.4660 0.6740 0.1880];
Rcolorcode=[0.4940 0.1840 0.5560];

SimMusclename=["bflh_r","bfsh_r","gaslat_r","gasmed_r","grac_r","recfem_r","sart_r","semimem_r","semiten_r","soleus_r","tfl_r","vasint_r","vaslat_r","vasmed_r"];
ExpMuscle=["RBICF","RSEMT","RMGAS","RRECF","RVASL","RVASM"];
SelectedMuscle=["RSEMT","RMGAS","RVASL","RRECF"];
%%% Defining transfer variable in the order to convert name of muscles in simulation to EMG sensors name 
%%% It must be in same dimention of ExpMuscle and each elemt of the Transfername is equal name to ExpMuscle
%%% Forexample ExpMuscle(1) is "RBICF" and equal name of this muscle in Rajagopal model would be "bflh_r" 
Transfername=["bflh_r","semiten_r","gasmed_r","recfem_r","vaslat_r","vasmed_r"];
RefT=0:.1:100;

%% Finding maximum EMG
SubjectNumber='T002';
Project='P005';
folder=append('C:\MyCloud\OneDriveUcf\Real\Simulation\Source\',SubjectNumber);
psname=append(Project,'_',SubjectNumber);
results_folder = append(folder,"\Result");   
load (append(results_folder,"\",psname,"_SimExpResultData.mat"));
%%

% xlabel(kk,'Normalized Time (%)');


if experflag
for A=1:length(AnalyzeMethod)
    for m=1:length(Modelname)
        results_folder = append(folder,"Result\",Modelname(m),"\");
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        for T1=1:length(Terials1)
            figure;
            movegui('center');
            kk = tiledlayout(4,4);
            title(kk,Terials1(T1));
            figure;
            movegui('center');
            t = tiledlayout(4,2);
            title(t,[Terials1(T1) 'Vecocity-Force-Motion']);
            TrialCounter=0;
            count=0;
            for T2=6:length(Terials2)
                
                TrialCounter=TrialCounter+1;
                filename=append(Terials1(T1),"_",Terials2(T2));
                % Reading Data
                ExpMotion=(ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).NormalMotion)./pi()*180; %scaling
                Time=ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).time.NormalT;
                ExpForce=ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).NormalExpForce;
                % Calculating Velocity 
                Vel=[0 0 0;diff(ExpMotion)./diff(Time)];
                % Making avarage and +- sd of Motion 
                RA1=mean(ExpMotion,2)-std(ExpMotion,0,2);
                RA2=mean(ExpMotion,2)+std(ExpMotion,0,2);
                RA3=mean(ExpMotion,2);
                % Making avarage and +- sd of Force              
                LA1=mean(ExpForce,2)-std(ExpForce,0,2);
                LA2=mean(ExpForce,2)+std(ExpForce,0,2);
                LA3=mean(ExpForce,2);
                %%% ploting motion 
                nexttile(t);
                x=RefT;
                x2 = [x, fliplr(x)];
                RABetween = [RA1', fliplr(RA2')];
                colororder([Lcolorcode;Rcolorcode])
                yyaxis left
                %%% ploting std of motion 
                fill(x2, RABetween,[Lcolorcode]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                hold on
                %%% ploting avarage of motion
                plot(x,RA3,'color',[Lcolorcode],'LineStyle','-');
                ylabel('Angle(degree)');
                %%% changing y axis
                yyaxis right
                LABetween = [LA1', fliplr(LA2')];
                %%% ploting std of force 
                fill(x2, LABetween,[Rcolorcode]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                %%% ploting avarage of force
                plot(x,LA3,'color',[Rcolorcode],'LineStyle','-');
                xlabel('Normalized Time (%)');
                ylabel('Force(N)');
                hold off
                %%% ploting Velocity
                nexttile;
                VA1=mean(Vel,2)-std(Vel,0,2);
                VA2=mean(Vel,2)+std(Vel,0,2);
                VA3=mean(Vel,2);
                VABetween = [VA1', fliplr(VA2')];
                fill(x2, VABetween,[Lcolorcode]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                hold on
                plot(x,VA3,'color',[Lcolorcode],'LineStyle','-');
                ylabel('Veocity(degree/s)');
                %%% Ploting EMG and activation of each muscle              
                for Flexmus=1:length(SelectedMuscle)
                    count=count+1;
                    % Making new plot
                    nexttile(kk)
                    colororder([Lcolorcode;Rcolorcode])
                    % preapring data
                    ExpMucs=ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).NormalExpEMG.(SelectedMuscle(Flexmus)); %scaling
%                     for i=1:3
%                     y2=interp1(ExpMotion(:,i),ExpMucs(:,i),x','linear','extrap');
%                     end
                    SimMusc=ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).NormalSimActivation.(Transfername((SelectedMuscle(Flexmus)==ExpMuscle)));
                    RA1=mean(ExpMucs,2)-std(ExpMucs,0,2);
                    RA2=mean(ExpMucs,2)+std(ExpMucs,0,2);
                    RA3=mean(ExpMucs,2);
                    LA1=mean(SimMusc,2)-std(SimMusc,0,2);
                    LA2=mean(SimMusc,2)+std(SimMusc,0,2);
                    LA3=mean(SimMusc,2);
                    x2 = [x, fliplr(x)];
                    RABetween = [RA1', fliplr(RA2')];
                    yyaxis right
                    fill(x2, RABetween,[Rcolorcode]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                    hold on
                    plot(x,RA3,'color',[Rcolorcode],'LineStyle','-');
                    ylabel('EMG');
                    % ylabel('Experiment activation');
                    hold on
                    yyaxis left
                    LABetween = [LA1', fliplr(LA2')];
                    fill(x2, LABetween,[Lcolorcode]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                    plot(x,LA3,'color',[Lcolorcode],'LineStyle','-');
                    if count<5
                        title([SelectedMuscle(Flexmus)]);
                    elseif count>12
                         xlabel('Normalized Time (%)');
                    elseif count
                    end
                    ylabel('activation ');
%                     legend('SD of Exp','SD of Sim','Mean of Exp','Mean of Sim');
                    hold off
                end
            end
        end
    end
end
end
% Trials=cell(length(Terials1).*length(Terials2),1);
%plotting avarage EMG 
 Trials={};
for T1=1:length(Terials1)
    for T2=1:length(Terials2)
        TrialCounter=TrialCounter+1;
        filename=append(Terials1(T1),"_",Terials2(T2));
        for Flmus=1:length(ExpMuscle)
             for itr=1:length(fieldnames(ResultData.(filename).time.Exp))
            
            MeanEmg(itr)=mean(ResultData.(filename).ExpEMG.(ExpMuscle(Flmus)).(Terials3(itr)));
            MeanActivation(itr)=mean(ResultData.(filename).(Modelname(1)).(AnalyzeMethod(1)).SimActivation.(Transfername(Flmus)).(Terials3(itr)));
             end
             MeanEmg2(TrialCounter,Flmus)=mean(MeanEmg,2);
             MeanActivation2(TrialCounter,Flmus)=mean(MeanActivation,2);
        end
        Trials=[Trials;filename];
    end
end
% T = array2table(MeanEmg2,'VariableNames',ExpMuscle);
% T.Trails=Trials;
% writetable(T,append(results_folder,'\EMG.csv'))
% T = T(:,ExpMuscle);
MaxEMG=max(MeanEmg2);

%%%%% EMG plot
 a=["Ex_IsoM","Fl_IsoM","Ex_IsoK","Fl_IsoK",];
 figure
 tid = tiledlayout(2,2);
for e=1:1:2
nexttile(e);
Y=MeanEmg2([1:5]+((e-1)*9),:)./MaxEMG;
X=[10 30 45 60 90];
plot(X,Y,'-*');
title(a(e));
ylabel('EMG');
ylim([0 1])
legend(ExpMuscle,'Location','northwest')
end
for e=1:1:2
nexttile(e+2);
Y=MeanEmg2([6:9]+((e-1)*9),:)./MaxEMG;
X=[60 120 180 240];
plot(X,Y,'-*');
title(a(e+2));
ylabel('EMG');
legend(ExpMuscle,'Location','northwest')
ylim([0 1])
end
title(tid,'Experiment')

%%%%%% Activation plot
a=["Ex_IsoM","Fl_IsoM","Ex_IsoK","Fl_IsoK",];
 figure
 tid = tiledlayout(2,2);
for e=1:1:2
nexttile(e);
Y=MeanActivation2([1:5]+((e-1)*9),:);
X=[10 30 45 60 90];
plot(X,Y,'-*');
title(a(e));
ylabel('Activation');
ylim([0 1])
legend(ExpMuscle,'Location','northwest')
end
for e=1:1:2
nexttile(e+2);
Y=MeanActivation2([6:9]+((e-1)*9),:);
X=[60 120 180 240];
plot(X,Y,'-*');
title(a(e+2));
ylabel('Activation');
ylim([0 1])
legend(ExpMuscle,'Location','northwest')
end
title(tid,'Simulation')

%%%%%%%%%
figure
tid = tiledlayout(1,2);
%%
a=["IsoM","IsoK"];
ratio=MeanEmg2([1:8],:)./MeanEmg2([9:16],:);
qq=0;
for kk=1:4:8
qq=qq+1;
nexttile(qq);
% Y=ResultData.EMG.data([kk:kk+3],:);
X=[10 30 45 60 90];
plot(X,ratio([kk:kk+3],:),'-*');
legend(ExpMuscle)
title(a(qq));
end
legend(ExpMuscle)
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

