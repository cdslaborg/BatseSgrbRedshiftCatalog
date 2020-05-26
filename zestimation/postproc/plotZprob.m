close all;
%clear all;
format compact; format long;
filePath = mfilename('fullpath');
addpath(genpath("../lib/"));
addpath(genpath('../../../../../lib/matlab/')) % lib codes

% change directory to the srouce code directory
filePath = mfilename('fullpath');
[scriptPath,~,~] = fileparts(filePath); cd(scriptPath); % Change working directory to source code directory.
cd(scriptPath); % Change working directory to source code directory.

figExportRequested = true;
fontSize = 12;

%rateModelList = ["H06","L08","B10","M14","M17","F18"];
rateModelList = "L08";
for rateModel = rateModelList

    %disp("processing rate model " + rateModel)

    buildDir = "D:\Dropbox\Projects\20181213_BatseLgrbRedshift\git\zestimation\build\winx64\intel\19.0.4.245\release\static\heap\serial\fortran\kfacOneThird\";
    rootDir = buildDir + rateModel + "\bin\out";

    zprobFileList = dir(fullfile(rootDir,'zprob_*'));
    zprobFileList = string({zprobFileList(:).name});
    zprobFileListLen = length(zprobFileList);

    z = importdata(fullfile(rootDir,'zgrid.txt'));

    figure; hold on; box on;
    for i = 1:25:zprobFileListLen
        d = importdata(fullfile(rootDir,zprobFileList(i)));
        plot( z.data(:,1) ...
            , d.data(:,1) ...
            ..., 'linewidth', 1 ...
            );
    end
    xlim([0 5]);
    xlabel("Redshift ( z )","interpreter","tex");
    ylabel("Probability Density Function");

    set(gca,'fontSize', fontSize)
    set(gca,'color', 'white')
    set(gcf,'color','white')
    if figExportRequested
        export_fig(fullfile(scriptPath, "zprob_" + rateModel + ".png"),"-m4 -transparent")
        hold off; close(gcf);
    else
        hold off;
    end

    %set(gca,'xscale','log','yscale','linear')

    zstat = importdata(fullfile(rootDir,'batse_zstat.txt'));
    interval50 = zstat.data(:,7) - zstat.data(:,5);
    interval90 = zstat.data(:,8) - zstat.data(:,4);
    disp( rateModel + ": " + string( sprintf('%0.2f',mean(interval50))) + ", " + string( sprintf('%0.2f',mean(interval90)) ) );
    %mean(zstat.data(:,7) ./ zstat.data(:,5))
    %mean(zstat.data(:,8) ./ zstat.data(:,4))

end

