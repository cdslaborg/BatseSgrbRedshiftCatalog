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

if exist('d','var')
    warning("skipping input data reading...");
else
    d = importdata('batseData.xlsx');
end

icolP64 = 19;
icolP1024 = 25;
icolDur = 7;

% find best fit parameters

logDur = log([ d.data.SGRBs(:,7) ; d.data.LGRBs(:,7) ]);
logRatioPF = log([ d.data.SGRBs(:,19) ./ d.data.SGRBs(:,25) ; d.data.LGRBs(:,19) ./ d.data.LGRBs(:,25) ]);

getFunc = @(param,logDur) param(1) + param(2) * erfc( (logDur-param(3)) / param(4) );
param0 = [0.077,0.24,-0.72,0.74];
%getFunc = @(param,logDur) param(1) + (log(16)/2) * erfc( (logDur-param(2)) / param(3) );
%param0 = [0.077,-0.72,0.74];
%getFunc = @(param,logDur) param(1) + (log(16)/2) * erfc( (logDur-log(0.256)) / param(2) );
%param0 = [0.077,0.74];

[paramBest,resnorm,~,exitflag,output] = lsqcurvefit(getFunc,param0,logDur,logRatioPF);

% plot results

figure; hold on; box on;

    dset  = ["LGRBs","SGRBs"];
    dsetColor = ["red","blue"];
    for i = 1:length(dset)
        plot( d.data.(dset(i))(:,7) ...
            , d.data.(dset(i))(:,19) ./ d.data.(dset(i))(:,25) ...
            , "." ...
            , "color", dsetColor(i) ...
            );
    end

    limDur = [0.005,2000];
    limRatioPF = [0.8,30];

    rangeLogDur = log(limDur(1)):0.01:log(limDur(2));
    rangeLogRatioPF = getFunc(paramBest,rangeLogDur);
    plot(exp(rangeLogDur), exp(rangeLogRatioPF),'linewidth',2,'color','black');

    set(gca,'xscale','log','yscale','log');
    xlim(limDur);
    ylim(limRatioPF);

hold off;

