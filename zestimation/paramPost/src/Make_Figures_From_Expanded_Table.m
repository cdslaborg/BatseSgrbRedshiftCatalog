% Created by Rebecca Proni on June 30th 2018
%
% Input files: AMHMCMC_DRBL_expanded.mat
% Output files: Some histograms and some plotted autocorrelation functions
%
% Plots the histograms and autocorrelation functions of the third to the
% twentieth variables contained in the table AMHMCMC_DRBL_expanded.mat, and
% saves these figures as png files

%clear all;
close all;
format compact; format long; fontSize = 13;
filePath = mfilename('fullpath');
[scriptPath,filePath,fileExt] = fileparts(filePath); cd(scriptPath);
%addpath(genpath('../../lib/matlab/')) % Amir's code lib path
%addpath(genpath('../../../MATLAB/library/')) % Rebecca's code lib path

Path.root = './';
Path.in   = '../tab/';
Path.out  = '../figures/';

%myTable = importdata([Path.in,thisModel,'.mat']);

Model.Type = {'AMHMCMC_DRBL_expanded'};
myTable = importdata([Path.in,thisModel,'.mat']);
thisModel = Model.Type;
thisModel = thisModel{1};

variableCellArray = {
                     myTable.mean_Log_Liso,  myTable.mean_Log_Epkz...
                     myTable.mean_Log_Eiso,  myTable.mean_Log_T90z...
                     myTable.Log_sigma_Liso, myTable.Log_sigma_Epkz...
                     myTable.Log_sigma_Eiso, myTable.Log_sigma_T90z...
                     myTable.rho_Liso_Epkz,  myTable.rho_Liso_Eiso...
                     myTable.rho_Liso_T90z,  myTable.rho_Eiso_Epkz...
                     myTable.rho_Epkz_T90z,  myTable.rho_Eiso_T90z...
                     myTable.mean_LogThresh, myTable.Log_ScaleThresh...
                    };
                
titleNamesCellArray = {
                       'Mean Log(L_{iso})',    'Mean Log(E_{pkz})'...
                       'Mean Log(E_{iso}',     'Mean Log(T_{90z})'...
                       'Log(\sigma L_{iso})',  'Log(\sigma E_{pkz})'...
                       'Log(\sigma E_{iso})',  'Log(\sigma T_{90z})'...
                       '\rho L_{iso} E_{pkz}', '\rho L_{iso} E_{iso}'...
                       '\rho L_{iso} T_{90z}', '\rho E_{iso} E_{pkz}'...
                       '\rho E_{pkz} T_{90z}', '\rho E_{iso} T_{90z}'...
                       'Mean Log(threshold)',  'Log(scale threshold)'
                      };
                  
fileNamesCellArray = {
                      'mean_log_liso',  'mean_log_Epkz'...
                      'mean_Log_Eiso',  'mean_Log_T90z'...
                      'Log_sigma_Liso', 'myTable.Log_sigma_Epkz'...
                      'Log_sigma_Eiso', 'Log_sigma_T90z'...
                      'rho_Liso_Epkz',  'rho_Liso_Eiso'...
                      'rho_Liso_T90z',  'rho_Eiso_Epkz'...
                      'rho_Epkz_T90z',  'rho_Eiso_T90z'...
                      'mean_LogThresh', 'Log_ScaleThresh'...
                      };
    
for i = 1:length(variableCellArray)
    varToPlotName = variableCellArray{i};
    titleName     = titleNamesCellArray{i};
    fileName      = fileNamesCellArray{i};

    % Histograms
    fprintf('Constructing the histogram for:\n')
    display(titleName)
    figure('visible','off','Color','none'); hold on; box on;
        histogram(varToPlotName,'EdgeAlpha',0,'FaceAlpha',1)
        title(strcat('Histogram of',{' '},titleName),'Interpreter', 'tex');
        %ylabel('Frequencey')
        %xlabel('Log of ', 'Interpreter', 'tex', 'fontSize', fontSize)
        set(gca,'color','none', 'fontSize', fontSize)
    export_fig ([Path.out,strcat(fileName,'_hist.png')],'-m2 -transparent');
    hold off; close(gcf);

    % ACF figures
    fprintf('Constructing the acf for:\n')
    display(titleName)
    figure('visible','off','Color','none'); hold on; box on;
        autocorr(varToPlotName,'NumLags',length(varToPlotName)-1)
        title(strcat('Autocorrelation of',{' '},titleName),'Interpreter', 'tex')
        %ylabel('Frequencey')
        %xlabel('Log of ', 'Interpreter', 'tex', 'fontSize', fontSize)
        set(gca,'color','none', 'fontSize', fontSize,'xscale','log')
    export_fig ([Path.out,strcat(fileName,'_acf.png')],'-m2 -transparent');
    hold off; close(gcf);
end

disp('...done!')