import org.opensim.modeling.*
path='C:\Program Files\OpenSim 4.1\Geometry';
 ModelVisualizer.addDirToGeometrySearchPaths(path);
%%%%%%%%% Scale ToolBox %%%%%
% ----------------------------------------------------------------------- %

% myModel = modeling.Model("gait2354_simbody.osim");               % Load a Model from file
% markerSetFile = "gait2354_MarkerSet.xml";                % Define the full path to the MarkerSet file
% newMarkers = modeling.MarkerSet(myModel, markerSetFile); % Construct a MarkerSet Object
% myModel.updateMarkerSet(newMarkers);                 

% Go to the folder in the subject's folder where .trc files are
trc_folder = 'C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody';

% specify where results will be printed.
results_folder = 'C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody\Out_put';

% Get and operate on the files
genericSetupPath='C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody';
genericSetupForScale='subject01_Setup_Scale.xml';
dir=fullfile(genericSetupPath,genericSetupForScale);
% scale = ScaleTool(dir);
% scale.run();

%%%%%%%%%%%%%%%%% IK %%%%%%%%%%%%%
% Pull in the modeling classes straight from the OpenSim distribution



% % Get the model
% [modelFile,modelFilePath,FilterIndex] = ...
%     uigetfile('*.osim','Pick the the model file to be used.');
% 
% % Load the model and initialize
% model = Model(fullfile(modelFilePath, modelFile));
% model.initSystem();
% 
% % Tell Tool to use the loaded model
% ikTool.setModel(model);
% 
% trialsForIK = dir(fullfile(trc_data_folder, '*.trc'));
% 
% nTrials = size(trialsForIK);
% 
% % Loop through the trials
% for trial= 1:nTrials
%     
%     % Get the name of the file for this trial
%     markerFile = trialsForIK(trial).name;
%     
%     % Create name of trial from .trc file name
%     name = regexprep(markerFile,'.trc','');
%     fullpath = fullfile(trc_data_folder, markerFile);
%     
%     % Get trc data to determine time range
%     markerData = MarkerData(fullpath);
%     
%     % Get initial and intial time 
%     initial_time = markerData.getStartFrameTime();
%     final_time = markerData.getLastFrameTime();
%     
%     % Setup the ikTool for this trial
%     ikTool.setName(name);
%     ikTool.setMarkerDataFileName(fullpath);
%     ikTool.setStartTime(initial_time);
%     ikTool.setEndTime(final_time);
%     ikTool.setOutputMotionFileName(fullfile(results_folder, [name '_ik.mot']));
%     
%     % Save the settings in a setup file
%     outfile = ['Setup_IK_' name '.xml'];
%     ikTool.print(fullfile(genericSetupPath, outfile));
%     
%     fprintf(['Performing IK on cycle # ' num2str(trial) '\n']);
%     % Run IK
%     ikTool.run();
% 
% end
% Get the model
model_folder = 'C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody';
Setupfile='subject01_Setup_CMC.xml';
Scaledir=fullfile(model_folder, Setupfile);
% Load the model and initialize
% model = Model(fullfile(model_folder, Setupfile));
% model.initSystem();
% InverseKinematicsTool.setModel(model)
% Tell Tool to use the loaded model
cmc = CMCTool(Scaledir);
cmc.setDesiredKinematicsFileName('C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody\OutputReference\subject0_walk1_ik.mot'); 
cmc.setExternalLoadsFileName('C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody\subject01_walk1_grf.xml')
% look at AbstractTool () to find more subclass for time range 
cmc.setStartTime(0.9);
cmc.setFinalTime(1.1);
cmc.setResultsDir('C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody\CMC');
cmc.run();
% AbstractTool ()

