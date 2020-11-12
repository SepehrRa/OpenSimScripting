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

name='P001_T001';

genericSetupForScale='subject01_Setup_Scale.xml';

% Rajagopal_Biodex_reaction_force_P001_T000_I0.mot
% Rajagopal_Biodex_kinematics_P001_T000_I0_ Rajagopal.mot
AnalyzeMethod=["SOP","CMC"];
Modelname={'Rajagopal'};
Terials1=["Fl","Ex"];
Terials2=["IsoM10","IsoM30","IsoM60","IsoM90","IsoK60","IsoK120","IsoK180","IsoK240"];
%% Model initiate %%
% model = Model(Scalemodel);
% state = model.initSystem();


%%  Setup the External Force File  %%

TimeT=zeros(4,2);
% for A=1:length(AnalyzeMethod)
    for m=1:length(Modelname)
        results_folder = append(folder,"Result\",Modelname(m),"\");
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        CMCSetup=append('CMC\',Modelname(m),'\CMC_Setup_I0.xml');
        IDSetup=append('CMC\',Modelname(m),'\ID_Setup.xml');
        for T1=1:length(Terials1)
            for T2=1:length(Terials2)
                filename=append(Terials1(T1),"_",Terials2(T2));
                results_folder2=append(results_folder,AnalyzeMethod,"\",filename,"\");
                status = mkdir(results_folder2(1));
                status = mkdir(results_folder2(2));
                ExLoad=ExternalLoads(append(results_folder,"P005_T001_ExForce_RLeg.xml"),true);
                ExForcefile=append(folder,"Data\P005_T001_Rknee_",filename,"_Torque.mot");
                FTable=importdata(ExForcefile);
%                 [r,c]=find(strncmp(FTable.textdata,'reaction_torque_z',17));
%                 [indx,c]=find(abs(FTable.data(:,10))>35);
%                 Stime=FTable.data(indx([1;find(diff(indx)>10)]),1);
%                 Etime=FTable.data(indx([(find(diff(indx)>10));end]),1);
                Stime=FTable.data(1,1);
                Etime=FTable.data(end,1);
                % Stime=2;
                % Etime=2.1;
%                 TimeT(k,:)=[Stime(1),Etime(1)];
                ExLoad.setDataFileName(ExForcefile);
                NewExForcefile=append(results_folder,"ID\",filename,"_ExForce_Setup.xml");
                ExLoad.print(NewExForcefile)
                IkFile=append(folder,"Data\P005_T001_Rknee_",filename,"_Motion.mot");
%                 IKtable=importdata(IkFile);
                %% ID %%%%
                idTool=	InverseDynamicsTool(append(results_folder,"P005_T001_ID_Setup_ref.xml"));
                idTool.setStartTime(Stime(1));
                idTool.setEndTime(Etime(1));
                idTool.setCoordinatesFileName(IkFile);
                idTool.setExternalLoadsFileName(NewExForcefile);
                idTool.setResultsDir(append(results_folder,"ID\"))
                idTool.setOutputGenForceFileName(append(filename,"_ID.sto"))
                idTool.print(append(results_folder,"ID\",filename,"_ID_Setup.xml"));
                idTool.run();
                %% SOP %%%%%
                analysis = AnalyzeTool(append(results_folder,"P005_T001_SOP_Setup_ref.xml"));
                analysis.setInitialTime(Stime(1));
                analysis.setFinalTime(Etime(1));
                analysis.setLowpassCutoffFrequency(6);
                analysis.setCoordinatesFileName(IkFile);
                analysis.setExternalLoadsFileName(NewExForcefile);
                analysis.setLoadModelAndInput(true);
                analysis.setResultsDir(append(results_folder2(1)));
                analysis.print(append(results_folder2(1),filename,"_",AnalyzeMethod(1),"_Setup.xml"))
                analysis.run();
                %% CMC %%%%%
                cmc = CMCTool(append(results_folder,"P005_T001_CMC_Setup_ref.xml"));
                cmc.setName(append(Modelname(m),'_',filename,'_CMC'))
                cmc.setDesiredKinematicsFileName(IkFile);
                cmc.setExternalLoadsFileName(NewExForcefile);
                cmc.setStartTime(Stime(1));
                cmc.setFinalTime(Etime(1));
                cmc.setResultsDir(append(results_folder2(2)));
                cmc.print(append(results_folder2(2),filename,"_",AnalyzeMethod(2),"_Setup.xml"));
                cmc.run();
                clear cmc ExLoad FTable idTool analysis
            end
        end
    end
% end

% look at AbstractTool () to find more subclass for time range 

