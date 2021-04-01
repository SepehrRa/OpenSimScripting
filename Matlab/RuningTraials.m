clear all
import org.opensim.modeling.*
path='C:\Program Files\OpenSim 4.1\Geometry';
ModelVisualizer.addDirToGeometrySearchPaths(path);
%% File address %%
folder = 'C:\MyCloud\OneDriveUcf\Real\Simulation\Source\T002\';
psname='P005_T002';
IKSteup='subject01_Setup_IK.xml';
DynamicMarkerFile='New_subject01_walk1.trc';
Resultdir='\Result\CMC\Rajagopal\';
ExForceSetup='P005_T001_ExForce.xml';

name='P001_T001';

genericSetupForScale='subject01_Setup_Scale.xml';

AnalyzeMethod=["SOP","CMC"];
Modelname=["Arnold","Rajagopal"];
Terials1=["Fl","Ex"];
Terials2=["IsoM10","IsoM30","IsoM45","IsoM60","IsoM90","IsoK60","IsoK120","IsoK180","IsoK240"];
Terials3=["iter1","iter2","iter3"];
load (append(folder,"Result\",psname,"_ResultData.mat"));
%% Runing simulation
TimeT=zeros(4,2);
    for m=1:length(Modelname)
        results_folder = append(folder,"Result\",Modelname(m),"\");
        status = mkdir(append(results_folder,"ID\"));
        %   IK_file=append(results_folder,name,'_ik.mot');
        %   NewExForcesetup=append(results_folder,'New_subject01_walk1_grf.xml');
        CMCSetup=append('CMC\',Modelname(m),'\CMC_Setup_I0.xml');
        IDSetup=append('CMC\',Modelname(m),'\ID_Setup.xml');
        for T1=1:length(Terials1)
            for T2=1:length(Terials2)
                filename=append(Terials1(T1),"_",Terials2(T2));
                ExLoad=ExternalLoads(append(results_folder,psname,"_ExForce_RLeg.xml"),true);
                ExForcefile=append(folder,"Data\",psname,"_Rknee_",filename,"_Torque.mot");
                FTable=importdata(ExForcefile);
                ExLoad.setDataFileName(ExForcefile);
                NewExForcefile=append(results_folder,"ID\",filename,"_ExForce_Setup.xml");
                ExLoad.print(NewExForcefile)
                IkFile=append(folder,"Data\",psname,"_Rknee_",filename,"_Motion.mot");
                %% ID %%%%
%                 idTool=	InverseDynamicsTool(append(results_folder,"P005_T001_ID_Setup_ref.xml"));
%                 idTool.setStartTime(FTable.data(1,1));
%                 idTool.setEndTime(FTable.data(end,1));
%                 idTool.setCoordinatesFileName(IkFile);
%                 idTool.setExternalLoadsFileName(NewExForcefile);
%                 idTool.setResultsDir(append(results_folder,"ID\"))
%                 idTool.setOutputGenForceFileName(append(filename,"_ID.sto"))
%                 idTool.print(append(results_folder,"ID\",filename,"_ID_Setup.xml"));
%                 idTool.run();
                for itr=1:length(fieldnames(ResultData.(filename).time.Exp))
                    results_folder2=append(results_folder,AnalyzeMethod,"\",filename,"\",Terials3(itr),"\");
                    status = mkdir(results_folder2(1));
                    status = mkdir(results_folder2(2));
                    Stime=ResultData.(filename).time.Exp.(Terials3(itr))(1);
                    Etime=ResultData.(filename).time.Exp.(Terials3(itr))(end);
                    %% SOP %%%%%
                    analysis = AnalyzeTool(append(results_folder,psname,"_SOP_Setup_ref.xml"));
                    analysis.setName(append(Modelname(m),'_',filename,'_',Terials3(itr)))
                    analysis.setInitialTime(Stime(1));
                    analysis.setFinalTime(Etime(1));
                    analysis.setLowpassCutoffFrequency(6);
                    analysis.setCoordinatesFileName(IkFile);
                    analysis.setExternalLoadsFileName(NewExForcefile);
                    analysis.setLoadModelAndInput(true);
                    analysis.setResultsDir(append(results_folder2(1)));
%                     analysis.print(append(results_folder2(1),filename,"_",AnalyzeMethod(1),"_Setup.xml"))
%                     analysis.run();
                    %% CMC %%%%%
                    cmc = CMCTool(append(results_folder,psname,"_CMC_Setup_ref.xml"));
                    cmc.setName(append(Modelname(m),'_',filename,'_',Terials3(itr),'_CMC'))
                    cmc.setDesiredKinematicsFileName(IkFile);
                    cmc.setExternalLoadsFileName(NewExForcefile);
                    cmc.setStartTime(Stime(1));
                    cmc.setFinalTime(Etime(1));
                    if strncmp(Terials2(T2),"IsoK",4)
                    cmc.setTimeWindow(0.005)
                    end
                    cmc.setResultsDir(append(results_folder2(2)));
                    cmc.print(append(results_folder2(2),filename,"_",AnalyzeMethod(2),"_Setup.xml"));
                    cmc.run();
                    clear cmc ExLoad FTable idTool analysis
                end
            end
        end
    end


% look at AbstractTool () to find more subclass for time range 

