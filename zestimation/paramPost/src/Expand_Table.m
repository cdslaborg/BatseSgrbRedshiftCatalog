% written by Rebecca Proni on June 29th 2018
%
% Input Files: AMHMCMC_DRBL_compact.xlsx
% Output Files: AMHMCMC_DRBL_expanded.mat
%
% Expands the table AMHMCMC_DRBL_compact.xlsx by duplicating each row by
% the given number of iterations.
% Then sorts this data so it is in chronological order and outputs the
% resulting table as a mat file in ../figures/

clear all;
close all;
format compact; format long;
filePath = mfilename('fullpath');
[scriptPath,filePath,fileExt] = fileparts(filePath); cd(scriptPath);
%addpath(genpath('../../lib/matlab/')) % Amir's code lib path
addpath(genpath('../../../MATLAB/library/')) % Rebecca's code lib path

Path.root = './';
Path.in = '../tab/';
%Path.in = '../in/';
Path.out = '../tab/';

Model.Type = {'AMHMCMC_DRBL_compact'};
thisModel  = Model.Type;
thisModel  = thisModel{1};

myData = importdata([Path.in,thisModel,'.xlsx']);
numberOfIterations{length(myData.data(:,2)),1} = [];

%calculates the number of times each row must be duplicated
for i = 1:length(myData.data(:,2))
    if i == length(myData.data(:,2))
        break
    end
    numberOfIterations{i} = myData.data(i+1,2) - myData.data(i,2);
end
numberOfIterations{end} = 0;

% Making the table
myTable = table();
myTable.variableCount   = myData.data(:,1);
myTable.numOfIterations = numberOfIterations;
myTable.mean_Log_Liso   = myData.data(:,5);
myTable.mean_Log_Epkz   = myData.data(:,6);
myTable.mean_Log_Eiso   = myData.data(:,7);
myTable.mean_Log_T90z   = myData.data(:,8);
myTable.Log_sigma_Liso  = myData.data(:,9);
myTable.Log_sigma_Epkz  = myData.data(:,10);
myTable.Log_sigma_Eiso  = myData.data(:,11);
myTable.Log_sigma_T90z  = myData.data(:,12);
myTable.rho_Liso_Epkz   = myData.data(:,13);
myTable.rho_Liso_Eiso   = myData.data(:,14);
myTable.rho_Liso_T90z   = myData.data(:,15);
myTable.rho_Eiso_Epkz   = myData.data(:,16);
myTable.rho_Epkz_T90z   = myData.data(:,17);
myTable.rho_Eiso_T90z   = myData.data(:,18);
myTable.mean_LogThresh  = myData.data(:,19);
myTable.Log_ScaleThresh = myData.data(:,20);

%{
% expanding the table- this part may take a few minutes

tic

for i = 1:length(myTable.variableCount)
    numOfTimesToDuplicate = myTable.numOfIterations{i}-1;
    for j = 1:numOfTimesToDuplicate
        myTable = [myTable;myTable(i,:)];
    end
    if mod((i+j),50) == 0
        fprintf('Have added the %1.fth row\n',i+j);
    end
end

toc

myTable = sortrows(myTable,1);
save('AMHMCMC_DRBL_expanded','myTable');

%}