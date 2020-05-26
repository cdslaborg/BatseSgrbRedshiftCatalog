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

if exist('belczynskiDelayedMergerRate','var')
    warning("skipping input data reading...");
else
    load('belczynskiDelayedMergerRate.mat');
end

% find best fit parameters

delay = belczynskiDelayedMergerRate(2:end-1,1);
count = exp( log(10) * belczynskiDelayedMergerRate(2:end-1,2) );

%getFunc = @(param,t) param(1) * lognpdf(t,0.1,param(2)); % lognormal
%param0 = [7.715888876935963   1.304336236656980];

mu = log(0.1);
getFunc = @(param,t) ( param(1) * lognpdf(t,log(0.1),param(2)) ); % lognormal
param0 = [5000 1.304336236656980];
[paramBest,resnorm,~,exitflag,output] = lsqcurvefit(getFunc,param0,delay,count);

% plot results

figure; hold on; box on;

    logDelayDif = max(log(delay)) - min(log(delay));
    logDelayLim = [ min(log(delay))-1*logDelayDif, max(log(delay))+1*logDelayDif ];
    logDelayRange = logDelayLim(1):0.01:logDelayLim(2);
    delayRange = exp(logDelayRange);

    plot( delayRange, getFunc(paramBest,delayRange),'linewidth',2,'color','black');

    plot( delay ...
        , count ...
        ..., "." ...
        ..., "markersize", 10 ...
        , 'linewidth', 2 ...
        , "color", "red" ...
        );

    set(gca,'xscale','log','yscale','log');
    %xlim(limDur);
    %ylim(limRatioPF);

hold off;

% resulting paramBest = [561.238487548257,0.961281252247549]
