close all;clc;
%clear all;
%clear classes;
format compact; format long;
%addpath(genpath('../../../../../libmatlab/')) % lib codes
addpath(genpath('../../../lib/paramonte/')) % ParaMonte lib codes




fontSize = 17;
out = 'DelayDistributionFigures';
outDir = [out,'/'];
exist(outDir,'dir')
if ~exist(outDir,'dir')
    outDir
    mkdir(outDir)
end

Model.ID = {'L08','H06', 'B10', 'M17', 'F18'}; %,
Model.count = length(Model.ID);


for i=1:Model.count
    figure()
    hold on
    StringModel=string(Model.ID(i));
    in='../../../mergerDelayDist/build/winx64/intel/19.1.1.216/release/static/serial/'+StringModel+'/romberg/bin/mergerDelayRate'+StringModel+'.txt';
    models=importdata(in);
    plot(models.data(:,1),models.data(:,2))
    plot(models.data(:,1),models.data(:,3),'linestyle','--')
    legend(StringModel+' delay',StringModel+' SFR')
    saveas(gcf,outDir+StringModel+'.png')
end