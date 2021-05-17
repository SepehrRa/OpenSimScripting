clear all
close all
% folder = 'C:\MyCloud\OneDriveUcf\Real\Simulation\Source\T002\';
AnalyzeMethod=["CMC"];
Modelname=["Arnold","Rajagopal"];
Terials1=["Fl"];
ExTerials1=["Fl","Ex"];
Terials2=["IsoM10","IsoM30","IsoM45","IsoM60","IsoK60","IsoK120"];
Terials3=["iter1","iter2","iter3"];
TimeT=zeros(4,2);
delimiterIn = ',';
M_ThresholdMin=20*3.14/180;
M_ThresholdMax=70*3.14/180;
ForceRatio=0.8;
TrialCounter=0;
readflage=1;
muscleplot=1;
experflag=1;
Lcolorcode=[0.4660 0.6740 0.1880];
Rcolorcode=[0.4940 0.1840 0.5560];
Pcolorcode=[1 0.2 0];

SimMusclename=["bflh_r","bfsh_r","gaslat_r","gasmed_r","recfem_r","sart_r","semimem_r","semiten_r","soleus_r","tfl_r","vasint_r","vaslat_r","vasmed_r"];
ExpMuscle=["RBICF","RSEMT","RMGAS","RRECF","RVASL","RVASM"];
%%% Defining transfer variable in the order to convert name of muscles in simulation to EMG sensors name 
%%% It must be in same dimention of ExpMuscle and each elemt of the Transfername is equal name to ExpMuscle
%%% Forexample ExpMuscle(1) is "RBICF" and equal name of this muscle in Rajagopal model would be "bflh_r" 
Transfername=["bflh_r","semiten_r","gasmed_r","recfem_r","vaslat_r","vasmed_r"];
%%% muscle for plot
SelectedMuscle=["RSEMT","RMGAS","RRECF","RVASL"];

RefT=0:.1:100;

SubjectNumber='T002';
Project='P005';
folder=append('C:\MyCloud\OneDriveUcf\Real\Simulation\Source\',SubjectNumber);
psname=append(Project,'_',SubjectNumber);
results_folder = append(folder,"\Result");   
load (append(results_folder,"\",psname,"_SimExpResultData.mat"));

%% Isokinetic
%%% Finding Max EMG
 Trials={};

 for T1=1:2
    for T2=1:6
        TrialCounter=TrialCounter+1;
        filename=append(ExTerials1(T1),"_",Terials2(T2));
        for Flmus=1:length(ExpMuscle)
             for itr=1:length(fieldnames(ResultData.(filename).time.Exp))
            
            MeanEmg(itr)=max(ResultData.(filename).ExpEMG.(ExpMuscle(Flmus)).(Terials3(itr)));
             end
             MeanEmg2(TrialCounter,Flmus)=mean(MeanEmg,2);
          
        end
        Trials=[Trials;filename];
    end
 end
 TrialCounter=0;
 for T1=1:1
    for T2=1:4
        TrialCounter=TrialCounter+1;
        filename=append(ExTerials1(T1),"_",Terials2(T2));
        for Flmus=1:length(ExpMuscle)
             for itr=1:length(fieldnames(ResultData.(filename).time.Exp))
            MeanArnoldActivation(itr)=mean(ResultData.(filename).(Modelname(1)).(AnalyzeMethod(1)).SimActivation.(Transfername(Flmus)).(Terials3(itr)));
            MeanRajActivation(itr)=mean(ResultData.(filename).(Modelname(2)).(AnalyzeMethod(1)).SimActivation.(Transfername(Flmus)).(Terials3(itr)));
             end
             MeanArnoldActivation2(TrialCounter,Flmus)=mean(MeanArnoldActivation,2);
             MeanRajActivation2(TrialCounter,Flmus)=mean(MeanRajActivation,2);
        end
        Trials=[Trials;filename];
    end
end
% T = array2table(MeanEmg2,'VariableNames',ExpMuscle);
% T.Trails=Trials;
% writetable(T,append(results_folder,'\EMG.csv'))
% T = T(:,ExpMuscle);

MaxEMG=max(MeanEmg2);

    if experflag
for A=1:length(AnalyzeMethod)
    for m=1
        results_folder = append(folder,"Result\",Modelname(m),"\");
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        for T1=1:length(Terials1)
            figure;
            movegui('center');
            kk = tiledlayout(2,2);
            title(kk,'Isokinetic trials');
            figure;
            movegui('center');
            t = tiledlayout(2,2);
            title(t,'Velocity-Force-Motion');
            TrialCounter=0;
            count=0;
            for T2=5:length(Terials2)
                
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
                Force1=mean(ExpForce,2)-std(ExpForce,0,2);
                Force2=mean(ExpForce,2)+std(ExpForce,0,2);
                Force3=mean(ExpForce,2);
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
                LABetween = [Force1', fliplr(Force2')];
                %%% ploting std of force 
                fill(x2, LABetween,[Rcolorcode]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                %%% ploting avarage of force
                plot(x,Force3,'color',[Rcolorcode],'LineStyle','-');
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
                
                for Flexmus=1:2
                    count=count+1;
                    % Making new plot
                    nexttile(kk)
                    colororder([Lcolorcode;Rcolorcode])
                    % preapring data
                    ExpMucs=ResultData.(filename).(Modelname(m)).(AnalyzeMethod(A)).NormalExpEMG.(SelectedMuscle(Flexmus))./MaxEMG((SelectedMuscle(Flexmus)==ExpMuscle)); %scaling
%                     for i=1:3
%                     y2=interp1(ExpMotion(:,i),ExpMucs(:,i),x','linear','extrap');
%                     end
                    SimMuscArnold=ResultData.(filename).(Modelname(1)).(AnalyzeMethod(A)).NormalSimActivation.(Transfername((SelectedMuscle(Flexmus)==ExpMuscle)));
                    SimMuscRaja=ResultData.(filename).(Modelname(2)).(AnalyzeMethod(A)).NormalSimActivation.(Transfername((SelectedMuscle(Flexmus)==ExpMuscle)));
                    EMG1=mean(ExpMucs,2)-std(ExpMucs,0,2);
                    EMG2=mean(ExpMucs,2)+std(ExpMucs,0,2);
                    EMG3=mean(ExpMucs,2);
                    Arnold1=mean(SimMuscArnold,2)-std(SimMuscArnold,0,2);
                    Arnold2=mean(SimMuscArnold,2)+std(SimMuscArnold,0,2);
                    Arnold3=mean(SimMuscArnold,2);
                    Raj1=mean(SimMuscRaja,2)-std(SimMuscRaja,0,2);
                    Raj2=mean(SimMuscRaja,2)+std(SimMuscRaja,0,2);
                    Raj3=mean(SimMuscRaja,2);
                    %%%
                    x2 = [x, fliplr(x)];                 
                    
                    hold on
                    plot(x,EMG3,'color',[Rcolorcode],'LineStyle','-','LineWidth',1);
                       plot(x,Arnold3,'color',[Lcolorcode],'LineStyle','-','LineWidth',1);
                    plot(x,Raj3,'color',[Pcolorcode],'LineStyle','-','LineWidth',1);
                    ylim([0 1])
                    %%%
                    EMGBetween = [EMG1', fliplr(EMG2')];
                    fill(x2, EMGBetween,[Rcolorcode]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                    ArnoldBetween = [Arnold1', fliplr(Arnold2')];
                    fill(x2, ArnoldBetween,[Lcolorcode]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                    RajBetween = [Raj1', fliplr(Raj2')];
                    fill(x2, RajBetween,[Pcolorcode]*.8,'EdgeAlpha',0,'FaceAlpha',0.2)
                    if count==1 || count==3
                        ylabel('Activation');
                    end
                    if count==2
                        legend('Exp','Arnold','Rajagopal');
                    end
%                     ylim([0 1])
                    title([SelectedMuscle(Flexmus)]);
                    %%%
                    if count>2
                         xlabel('Normalized Time (%)');
                    else
                    end
                    hold off
                end
            end
        end
    end
end
end
TrialCounter=0;
% Normalizing EMG base on the maximum value of EMG : 
% 1- Averaging data within the each trials. 2- Getting maximum of all Isometric trials   
%% Isometric condition
for T1=1:2
    for T2=1:6
        TrialCounter=TrialCounter+1;
        filename=append(ExTerials1(T1),"_",Terials2(T2));
        for Flmus=1:length(ExpMuscle)
             for itr=1:length(fieldnames(ResultData.(filename).time.Exp))            
             MeanEmgIsoM(itr)=mean(ResultData.(filename).ExpEMG.(ExpMuscle(Flmus)).(Terials3(itr)));
             end
             MeanEmgIsoM2(TrialCounter,Flmus)=mean(MeanEmgIsoM,2);
          
        end
    end
end
MaxEMGIsoM=max(MeanEmgIsoM2);
figure
tid = tiledlayout(2,2);
title(tid,'Isometric trials');
for Flexmus=1:length(SelectedMuscle)
nexttile
Y=MeanEmgIsoM2(1:4,(SelectedMuscle(Flexmus)==ExpMuscle))./MaxEMGIsoM((SelectedMuscle(Flexmus)==ExpMuscle));
X=[10 30 45 60];

Y2=MeanArnoldActivation2(:,(SelectedMuscle(Flexmus)==ExpMuscle));
Y3=MeanRajActivation2(:,(SelectedMuscle(Flexmus)==ExpMuscle));

plot(X,Y','o','color',Rcolorcode);
hold on
plot(X,Y2','*','color',Lcolorcode);
plot(X,Y3','d','color',Pcolorcode);
hold off
title(ExpMuscle((SelectedMuscle(Flexmus)==ExpMuscle)));
ylim([0 1])
if Flexmus==2
legend('Exp','Arnold','Rajagopal')
end
if Flexmus>2
xlabel('Knee angle (deg)');
ylim([0 0.1])
end
if Flexmus==1||Flexmus==3
ylabel('Activation');
end
end

% %%%%% EMG plot
%  a=["Ex_IsoM","Fl_IsoM","Ex_IsoK","Fl_IsoK",];
%  figure
%  tid = tiledlayout(2,2);
% for e=1:1:2
% nexttile(e);
% Y=MeanEmg2([1:5]+((e-1)*9),:)./MaxEMG;
% X=[10 30 45 60 90];
% plot(X,Y,'-*');
% title(a(e));
% ylabel('EMG');
% ylim([0 1])
% legend(ExpMuscle,'Location','northwest')
% end
% for e=1:1:2
% nexttile(e+2);
% Y=MeanEmg2([6:9]+((e-1)*9),:)./MaxEMG;
% X=[60 120 180 240];
% plot(X,Y,'-*');
% title(a(e+2));
% ylabel('EMG');
% legend(ExpMuscle,'Location','northwest')
% ylim([0 1])
% end
% title(tid,'Experiment')
% 
% %%%%%% Activation plot
% a=["Ex_IsoM","Fl_IsoM","Ex_IsoK","Fl_IsoK",];
%  figure
%  tid = tiledlayout(2,2);
% for e=1:1:2
% nexttile(e);
% Y=MeanArnoldActivation2([1:5]+((e-1)*9),:);
% X=[10 30 45 60 90];
% plot(X,Y,'-*');
% title(a(e));
% ylabel('Activation');
% ylim([0 1])
% legend(ExpMuscle,'Location','northwest')
% end
% for e=1:1:2
% nexttile(e+2);
% Y=MeanArnoldActivation2([6:9]+((e-1)*9),:);
% X=[60 120 180 240];
% plot(X,Y,'-*');
% title(a(e+2));
% ylabel('Activation');
% ylim([0 1])
% legend(ExpMuscle,'Location','northwest')
% end
% title(tid,'Simulation')
% 
% %%%%%%%%%
% figure
% tid = tiledlayout(1,2);
% %%
% a=["IsoM","IsoK"];
% ratio=MeanEmg2([1:8],:)./MeanEmg2([9:16],:);
% qq=0;
% for kk=1:4:8
% qq=qq+1;
% nexttile(qq);
% % Y=ResultData.EMG.data([kk:kk+3],:);
% X=[10 30 45 60 90];
% plot(X,ratio([kk:kk+3],:),'-*');
% legend(ExpMuscle)
% title(a(qq));
% end
% legend(ExpMuscle)
% %%
% for kk=1:4:16
% qq=qq+1;
% nexttile(e);
% Y=ResultData.MaxActivation.data([kk:kk+3],[2,10,5,6,14,15]);
% X=[10 30 60 90];
% plot(X,Y,'-*');
% legend(ExpMuscle)
% title(a(qq));
% end
% 
% X=[10 30 60 90];
%  a=["Fl_IsoK","Fl_IsoM","Ex_IsoK","Ex_IsoM"'];
% e = tiledlayout(2,2);
%  qq=0;
% for kk=1:4:16
% qq=qq+1;
% nexttile(e);
% Y=ResultData.MaxActuateForce.data([kk:kk+3],[2,10,5,6,14,15]);
% bar(X,Y','stacked')
% legend(ExpMuscle)
% title(a(qq));
% end

