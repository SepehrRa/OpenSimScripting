clear all
import org.opensim.modeling.*
path='C:\Program Files\OpenSim 4.1\Geometry';
ModelVisualizer.addDirToGeometrySearchPaths(path);
% C:\MyCloud\OneDriveUcf\Real\JD\GasContribution\OpenSimModel\SimulationDataAndSetupFiles-4.0\Subject1UnAffectedleg.osim
%% File address %%
folder = 'C:\MyCloud\OneDriveUcf\Real\Simulation\Source\T001\';
% Scalemodel={'C:\MyCloud\OneDriveUcf\Real\JD\GasContribution\OpenSimModel\SimulationDataAndSetupFiles-4.0\Subject1UnAffectedleg.osim',...
%     'C:\MyCloud\OneDriveUcf\Real\JD\GasContribution\OpenSimModel\PatellaModel\Scaled_UnAffected.osim'};

IKSteup='subject01_Setup_IK.xml';
DynamicMarkerFile='New_subject01_walk1.trc';
Resultdir='\Result\CMC\Rajagopal\';
ExForceSetup='P005_T001_ExForce.xml';
% Markerfile=[folder 'CMC\Rajagopal'];
% Modelname={'Thelen';'Rajagopal'};

name='subject01';

genericSetupForScale='subject01_Setup_Scale.xml';

% Rajagopal_Biodex_reaction_force_P001_T000_I0.mot
% Rajagopal_Biodex_kinematics_P001_T000_I0_ Rajagopal.mot
AnalyzeMethod={'SOP'};
Modelname={'Rajagopal'};
Terials1=["Ex","Fl"];
Terials2=["IsoK60","IsoK120","IsoK180","IsoK240","IsoM10","IsoM30","IsoM60","IsoM90"];
%% Model initiate %%
% model = Model(Scalemodel);
% state = model.initSystem();


%%  Setup the External Force File  %%

TimeT=zeros(4,2);
for A=1:length(AnalyzeMethod)
    for m=1:length(Modelname)
        results_folder = append(folder,"\Result\",Modelname(m),"\");
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        CMCSetup=append('CMC\',Modelname(m),'\CMC_Setup_I0.xml');
        IDSetup=append('CMC\',Modelname(m),'\ID_Setup.xml');
        if Modelname(m)=="Rajagopal"
            Sw=2;
            Extention="_Rajagopal";
        else
            Sw=1;
            Extention="";
        end
        for T1=1:length(Terials1)
            for T2=1:length(Terials2)
                ExLoad=ExternalLoads(append(results_folder,"P005_T001_ExForce_RLeg.xml"),true);
                ExForcefile=append(folder,"Data\P005_T001_Rknee_",Terials1(T1),"_",Terials2(T2),"_Torque.mot");
                FTable=importdata(ExForcefile);
%                 [r,c]=find(strncmp(FTable.textdata,'reaction_torque_z',17));
                [indx,c]=find(abs(FTable.data(:,10))>30);
                Stime=FTable.data(indx([1;find(diff(indx)>10)]),1);
                Etime=FTable.data(indx([(find(diff(indx)>10));end]),1);
%                 Stime=FTable.data(indx(1),1);
%                 Etime=FTable.data(indx(end),1);
                % Stime=2;
                % Etime=2.1;
%                 TimeT(k,:)=[Stime(1),Etime(1)];
                ExLoad.setDataFileName(ExForcefile);
                NewExForcefile=append(results_folder,AnalyzeMethod(A),"\",Terials1(T1),Terials1(T2),"\ExForce_.xml");
                ExLoad.print(NewExForcefile)
                IkFile=append(folder,"Rajagopal_Biodex_kinematics_P001_T000_",Terials2(T2),Extention,".mot");
                IKtable=importdata(IkFile);
                %% ID %%%%
                idTool=	InverseDynamicsTool(append(folder,IDSetup));
                idTool.setStartTime(Stime(1));
                idTool.setEndTime(Etime(1));
                idTool.setCoordinatesFileName(IkFile);
                idTool.setExternalLoadsFileName(NewExForcefile);
                idTool.setResultsDir(append(results_folder,Terials(k)))
                idTool.setOutputGenForceFileName(append("ID_",Modelname(m),"_",Terials(k),".sto"))
                idTool.print(append(results_folder,Terials(k),"\ID_Setup_",Terials(k),".xml"));
                idTool.run();
                if AnalyzeMethod(A)=="SOP"
                    %% SOP %%%%%
                    %                 motion = Storage(IK_file);
                    %                 static_optimization = StaticOptimization(append(folder,CMCSetup));
                    %                 static_optimization.setStartTime(Stime);
                    %                 static_optimization.setEndTime(Etime);
                    %                 static_optimization.setUseModelForceSet(true);
                    %                 static_optimization.setUseMusclePhysiology(true);
                    %                 static_optimization.setActivationExponent(2);
                    %                 static_optimization.setConvergenceCriterion(0.0001);
                    %                 static_optimization.setMaxIterations(100)
                    model.addAnalysis(append(folder,CMCSetup))
                    analysis = AnalyzeTool(model);
                    analysis.setName(name);
                    analysis.setModel(model);
                    analysis.setInitialTime(motion.getFirstTime());
                    analysis.setFinalTime(motion.getLastTime());
                    analysis.setLowpassCutoffFrequency(6);
                    analysis.setCoordinatesFileName(IK_file);
                    analysis.setExternalLoadsFileName(NewExForcesetup);
                    analysis.setLoadModelAndInput(true);
                    analysis.setResultsDir([results_folder 'SOP']);
                    analysis.run();
                else
                    %% CMC %%%%%
                    cmc = CMCTool(append(folder,CMCSetup));
                    cmc.setName(append('CMC_',Modelname(m),'_',Terials(k)))
                    cmc.setDesiredKinematicsFileName(IkFile);
                    cmc.setExternalLoadsFileName(NewExForcefile);
                    cmc.setStartTime(Stime);
                    cmc.setFinalTime(Etime);
                    cmc.setResultsDir(append(results_folder,Terials(k)));
                    cmc.print(append(results_folder,Terials(k),"\CMC_Setup_",Terials(k),".xml"));
                    cmc.run();
                    clear cmc ExLoad FTable idTool ExLoad
                end
            end
        end
    end
end

% look at AbstractTool () to find more subclass for time range 

