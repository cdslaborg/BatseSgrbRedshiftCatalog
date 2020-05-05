clear all;
close all;
format compact; format long; fontSize = 13;
filePath = mfilename('fullpath');
[scriptPath,filePath,fileExt] = fileparts(filePath); cd(scriptPath);
%addpath(genpath('../../../MATLAB/library/'))
addpath(genpath('/Users/bproni/Documents/MATLAB/library/'))

Input.Path.root = './';
Input.Path.in   = '/Users/bproni/Documents/git/BatseSgrbRedshiftCatalog/cosmicRate/in/'; 
Input.Path.out  = '../out/';

%/Users/bproni/Documents/Research/Make_Figures/data_565.out
myFilePath_LGRB = [Input.Path.in,'BATSE_1366_LGRB_Pbol_Epk_Sbol(0.001,20000).txt'];
myFilePath_SGRB = [Input.Path.in,'BATSE_565_SGRB_bol64(0.001,20000)_data_6904.txt'];


fid       = fopen(myFilePath_LGRB);
LGRB_data = textscan(fid,'%f %f %f %f %f %f %f %f','HeaderLines',1,'CollectOutput',true);
LGRB_data = LGRB_data{:};
fid       = fclose(fid);

LGRB_data_table.log_Pbol = LGRB_data(:,2);
LGRB_data_table.log_Sbol = LGRB_data(:,3);
LGRB_data_table.log_Epk = LGRB_data(:,4);
LGRB_data_table.log_T90 = LGRB_data(:,8);

variableCellArray_LGRB = {
                          LGRB_data_table.log_Pbol, LGRB_data_table.log_Sbol...
                          LGRB_data_table.log_Epk, LGRB_data_table.log_T90...
                         };

fid       = fopen(myFilePath_SGRB);
SGRB_data = textscan(fid,'%f %f %f %f %f %f %f %f %f','HeaderLines',1,'CollectOutput',true);
SGRB_data = SGRB_data{:};
fid       = fclose(fid);

SGRB_data_table.log_Pbol = SGRB_data(:,6);
SGRB_data_table.log_Sbol = SGRB_data(:,7);
SGRB_data_table.log_Epk = SGRB_data(:,8);
SGRB_data_table.log_T90 = SGRB_data(:,9);

variableCellArray_SGRB = {
                          SGRB_data_table.log_Pbol, SGRB_data_table.log_Sbol...
                          SGRB_data_table.log_Epk, SGRB_data_table.log_T90...
                         };
                
titleNamesCellArray = {
                       'Bolometric Peak Flux: Log( P_{bol} [erg s^{-1} cm^{-2}] )', 'Bolometric Fluence: Log( S_{bol} [erg cm^{-2}] )'...
                       'Observed Peak Energy: Log( E_{pk} [keV] )', 'Observed Duration: Log( T_{90} [s] )'...
                      };
                  
fileNamesCellArray = {
                      'Hist_Log_Pbol',  'Hist_Log_Sbol'...
                      'Hist_Log_Epk',  'Hist_Log_T90'...
                      };

for i = 1:length(variableCellArray_SGRB)
    varToPlotName_SGRB = variableCellArray_SGRB{i};
    varToPlotName_LGRB = variableCellArray_LGRB{i};
    titleName = titleNamesCellArray{i};
    fileName = fileNamesCellArray{i};
    
    fprintf('Constructing the histogram for:\n')
    display(titleName)
    figure('visible','off','Color','none'); hold on; box on;
        h1 = histogram(varToPlotName_SGRB,'EdgeAlpha',1,'FaceAlpha',0.6);
        h2 = histogram(varToPlotName_LGRB,'EdgeAlpha',1,'FaceAlpha',0.6);
        h1.EdgeColor = "none";  h2.EdgeColor = "none";
        %h1.FaceAlpha = 0; h2.FaceAlpha = 0;
        h1.BinWidth = 0.2; h2.BinWidth = 0.2;
        ylabel('Counts')
        xlabel(titleName, 'Interpreter', 'tex', 'fontSize', fontSize)
        set(gca,'color','none', 'fontSize', fontSize)
        legend('SGRB','LGRB')
    export_fig ([Input.Path.out,fileName,'.png'],'-m2 -transparent');
    hold off; close(gcf);
end

disp('...done!')