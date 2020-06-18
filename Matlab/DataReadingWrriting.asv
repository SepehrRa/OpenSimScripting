% Add btk library to MATLAB path  -> https://code.google.com/archive/p/b-tk/downloads
% Add Matlab Opensim Pipeline Tools  to MATLAB path -> https://simtk.org/frs/?group_id=660
clear all
folder = 'C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody\';
fname = 'Static.c3d';
%% C3D file reading 
data = btk_loadc3d([folder fname]);

%% Generate .Trc
markerset={'LASI','RASI','LPSI','RPSI','LKNE','LTHI','LANK','LTIB','LTOE','LHEE','RKNE','RTHI','RANK','RTIB','RTOE','RHEE'};
DataRate=148;
% % To convert C3D of REAL LAb axis to Opensim Axis 
for s=1:3:lenght(markerset)*3
   markerpos(:,[s s+1 s+2])=[NewData(:,s+1),-1*NewData(:,s+2),NewData(:,s)]; 
end
generateTrcFile([name '_10_20_Marker.trc'],DataRate, markerpos, markerset);
