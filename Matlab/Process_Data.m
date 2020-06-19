% Add btk library to MATLAB path  -> https://code.google.com/archive/p/b-tk/downloads
% Add Matlab Opensim Pipeline Tools  to MATLAB path -> https://simtk.org/frs/?group_id=660
clear all
folder = 'C:\MyCloud\GitHub\OpenSimScripting\Matlab\Data\Gait2354_Simbody\';
fname = 'Static.c3d';
q=1;
%% C3D file reading 
data = c3d_getdata([folder fname]);
%% Generate .Trc for Marker set
Markerset=fieldnames(data.marker_data.Markers);
%%% Or any new lable. So you can change your lable based on your model. %%%
% Newmarkerlable={'LASI','RASI','LPSI','RPSI','LKNE','LTHI','LANK','LTIB','LTOE','LHEE','RKNE','RTHI','RANK','RTIB','RTOE','RHEE'};
%%% To convert C3D Vicon of REAL LAb axis to Opensim Axis 
MarkerData=data.marker_data.Time;
for i = 1:length(Markerset)
   MarkerData = [MarkerData data.marker_data.Markers.(Markerset{i})(:,2) data.marker_data.Markers.(Markerset{i})(:,3) data.marker_data.Markers.(Markerset{i})(:,1)];
end
generate_Marker_Trc(Markerset,MarkerData,data.marker_data.Info);

