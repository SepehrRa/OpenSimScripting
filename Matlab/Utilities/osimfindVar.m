

%%
clear; clc; close all

% gem
% activation

fileName = 'Test2.mot';
% fileName = 'subject01_walk1_controls.sto';
varData = importdata(fileName);
time = varData.data(:, 1);
Varname = 'activation';
Mycellarray=strfind(varData.colheaders,Varname);
myIndex = find(cellfun('size', Mycellarray, 1) > 0);
DesiredData = varData.data(:,myIndex);
savefile = [fileName,'.mat'];
plot(time,DesiredData);
hleg1=legend(varData.colheaders{myIndex});
set(hleg1,'Interpreter','none')
save(savefile, 'time', 'DesiredData');

% dataSize = size(varData.data);
% for headerNum = 1:dataSize(2)
% 
% if strfind(varData.colheaders{headerNum},Varname);
%     a=a+1;
%     myIndex(a) = headerNum;
% end
% end
% myIndex=find(ismember(varData.colheaders,Varname));

