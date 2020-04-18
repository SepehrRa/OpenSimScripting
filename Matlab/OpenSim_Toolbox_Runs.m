clear all
import org.opensim.modeling.*
path='C:\Program Files\OpenSim 4.1\Geometry';
ModelVisualizer.addDirToGeometrySearchPaths(path);

%% File address %%
folder = 'C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody\';
Scalemodel='subject01_scaledOnly.osim';
CMCSetup='subject01_Setup_CMC.xml';
IKSteup='subject01_Setup_IK.xml';
DynamicMarkerFile='New_subject01_walk1.trc';
ExForceSetup='subject01_walk1_grf.xml';
Markerfile='New_subject01_walk1.trc';
results_folder = [folder 'Out_put\'];
genericSetupForScale='subject01_Setup_Scale.xml';
name='subject01';

%% Scale ToolBox %%
scale = ScaleTool([folder genericSetupForScale]);
scale.run();

%% Model initiate %%
model = Model(fullfile([folder, Scalemodel]));
model.initSystem();

%%  Setup the External Force File  %%
ExLoad=ExternalLoads([folder ExForceSetup],true);
ExLoad.setDataFileName([folder 'New_subject01_walk1_grf.mot']);
ExLoad.print([results_folder 'New_subject01_walk1_grf.xml'])

%% IK %%
%%% Get trc data to determine time range
markerData = MarkerData([folder Markerfile]); 
initial_time = markerData.getStartFrameTime();
final_time = markerData.getLastFrameTime();
%%%
%%% To creat ik tool without any xml file
% ikTool=InverseKinematicsTool();
% ikTool.setName(name);
%%% 
%%% To creat ik tool with xml file and further changes
ikTool=InverseKinematicsTool([folder IKSteup]); % to read xml file for IK
ikTool.setModel(model);
ikTool.setMarkerDataFileName([folder DynamicMarkerFile]);
ikTool.setStartTime(initial_time);
ikTool.setEndTime(final_time);
ikTool.setOutputMotionFileName([results_folder name '_ik.mot']);
ikTool.print([results_folder name '_IK_Setup.xml']);
% ikTool.getIKTaskSet() % to change the scale files
ikTool.run();

%%  SOP %%

%%  CMC  %%

% look at AbstractTool () to find more subclass for time range 
cmc = CMCTool([folder CMCSetup]);
cmc.setDesiredKinematicsFileName([results_folder name '_ik.mot']);
cmc.setExternalLoadsFileName([results_folder 'New_subject01_walk1_grf.xml']);
cmc.setStartTime(0.9);
cmc.setFinalTime(0.95);
cmc.setResultsDir([folder 'CMC']);
cmc.run();
