% comparing old rival estimates with each other 
%close all;
%clear all;
format compact; format long;
addpath(genpath("../")) % local lib codes
addpath(genpath("../../../../../libmatlab/")) % lib codes

% change directory to the srouce code directory
filePath = mfilename("fullpath");
[scriptPath,~,~] = fileparts(filePath); cd(scriptPath); % Change working directory to source code directory.
cd(scriptPath); % Change working directory to source code directory.

inputDataRootPath = "D:\Dropbox\Projects\20181213_BatseSgrbRedshift\git\mergerDelayDist\build\winx64\intel\19.0.4.245\release\static\serial\B10\romberg\bin";
if exist('delayedMergerRate','var')
    warning("skipping input data reading...");
else
    delayedMergerRate = importdata(fullfile(inputDataRootPath,"mergerDelayRateB10.txt"));
    figure
    plot(delayedMergerRate.data(:,1),delayedMergerRate.data(:,3))
    hold on
    plot(delayedMergerRate.data(:,1),delayedMergerRate.data(:,2))
end

% polynomial fit to this function

global a
a = 100;
