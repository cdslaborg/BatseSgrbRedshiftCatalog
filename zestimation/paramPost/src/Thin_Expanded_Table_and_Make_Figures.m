% Created by Rebecca Proni on June 30th 2018
%
% Input files: myThinTable.mat
% Output files: Some histograms and some plotted autocorrelation functions
%
% Thins the tableby a given valueval
% Plots the autocorrelation functions of the third to the twentieth 
% variables contained in the table myThinTable.mat, and saves these figures
% as png files

clear all;
close all;
format compact; format long; fontSize = 13;
filePath = mfilename('fullpath');
[scriptPath,filePath,fileExt] = fileparts(filePath); cd(scriptPath);

Path.root = './';
Path.in   = '../tab/';
Path.out  = '../figures/';

Model.Type      = {'AMHMCMC_DRBL_expanded'};
thisModel       = Model.Type;
thisModel       = thisModel{1};
myImportedTable = importdata([Path.in,thisModel,'.mat']);

myNewTable    = table();
thinningValue = 150;
n = 1;

%thinning the table
for i = 1:floor(height(myImportedTable)/thinningValue)
    myNewTable = [myNewTable;myImportedTable(n,:)];
    n = n + thinningValue;
end

save([Path.in,'AMHMCMC_DRBL_thinned.mat'],'myNewTable');

variableCellArray = {
                     myNewTable.mean_Log_Liso,  myNewTable.mean_Log_Epkz...
                     myNewTable.mean_Log_Eiso,  myNewTable.mean_Log_T90z...
                     myNewTable.Log_sigma_Liso, myNewTable.Log_sigma_Epkz...
                     myNewTable.Log_sigma_Eiso, myNewTable.Log_sigma_T90z...
                     myNewTable.rho_Liso_Epkz,  myNewTable.rho_Liso_Eiso...
                     myNewTable.rho_Liso_T90z,  myNewTable.rho_Eiso_Epkz...
                     myNewTable.rho_Epkz_T90z,  myNewTable.rho_Eiso_T90z...
                     myNewTable.mean_LogThresh, myNewTable.Log_ScaleThresh...
                    };
                
titleNamesCellArray = {
                       'Mean Log(L_{iso})',    'Mean Log(E_{pkz})'...
                       'Mean Log(E_{iso})',     'Mean Log(T_{90z})'...
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
                      'Log_sigma_Liso', 'Log_sigma_Epkz'...
                      'Log_sigma_Eiso', 'Log_sigma_T90z'...
                      'rho_Liso_Epkz',  'rho_Liso_Eiso'...
                      'rho_Liso_T90z',  'rho_Eiso_Epkz'...
                      'rho_Epkz_T90z',  'rho_Eiso_T90z'...
                      'mean_LogThresh', 'Log_ScaleThresh'...
                      };
                      
% Plotting the histograms and ACF figures for thin tables
for i = 1:length(variableCellArray)
    varToPlotName = variableCellArray{i};
    titleName = titleNamesCellArray{i};
    fileName = fileNamesCellArray{i};
    
    fprintf('Constructing the histogram for:\n')
    display(titleName)
    figure('visible','off','Color','none'); hold on; box on;
        histogram(varToPlotName,'EdgeAlpha',0,'FaceAlpha',1)
        title(strcat('Histogram of',{' '},titleName,{' '},'thinned'),'Interpreter', 'tex');
        %ylabel('Frequencey')
        %xlabel('Log of ', 'Interpreter', 'tex', 'fontSize', fontSize)
        set(gca,'color','none', 'fontSize', fontSize)
    export_fig ([Path.out,strcat(fileName,'_hist_thinned.png')],'-m2 -transparent');
    hold off; close(gcf);

    fprintf('Constructing the acf for:\n')
    display(titleName)
    figure('visible','off','Color','none'); hold on; box on;
        autocorr(varToPlotName,'NumLags',length(varToPlotName)-1)
        title(strcat('Autocorrelation of',{' '},titleName,{' '},'thinned'),'Interpreter', 'tex')
        %ylabel('Frequencey')
        %xlabel('Log of ', 'Interpreter', 'tex', 'fontSize', fontSize)
        set(gca,'color','none', 'fontSize', fontSize,'xscale','log')
    export_fig ([Path.out,strcat(fileName,'_acf_thinned.png')],'-m2 -transparent');
    hold off; close(gcf);
end

disp('...done!')