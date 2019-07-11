%% osimfindVar finds desierd variable from .mot/.sto files
% Input  = .mot or sto files.
% Output = mat file of desired variables and ploting them.
%----------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and           %
% simulation. See http://opensim.stanford.edu and the NOTICE file         %
% for more information.
% This function finds any information about specific variable in .mot/.sto 
% fiels. First change the directory of Current Folder to the input files dir.
% By Sepehr Ramezani (2019)
%
%

%%
clear; clc; close all
%% Replace 'Test2.mot' with the file which you want to search in.
fileName = 'Test2.mot';
%% Replace 'gem' with the variable name which you want explore in the file.
Varname = 'gem';
varData = importdata(fileName);
time = varData.data(:, 1);
% searching variable
Mycellarray=strfind(varData.colheaders,Varname);
myIndex = find(cellfun('size', Mycellarray, 1) > 0);
DesiredData = varData.data(:,myIndex);
savefile = [fileName,'.mat'];
plot(time,DesiredData);
hleg1=legend(varData.colheaders{myIndex});
set(hleg1,'Interpreter','none')
save(savefile, 'time', 'DesiredData');

