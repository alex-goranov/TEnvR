%% Description

% This code is used to process raw EEM spectra that have been acquired on a Horiba Aqualog. 

% The spectra must have been pre-processed immediately after acquisition:
    % Instrument-specific correction (automatic)
    % Rayleigh masking 
    % Inner-filter effect correction
    % Normalization to Raman Units (RU) or Quinine Sulfate Units (QSU)

% Spectra for this version of the code must have been exported as .dat files.
% Note that the "readineems" code of drEEM is capable of importing EEM 
% files in various formats: 'csv', 'xls', 'xlsx','dat', 'txt' - modify
% below if needed. 

% The EEM_Process_Aqualog script is based on the tutorials by Murphy et al. 2013. 
% This script utilizes functions from the drEEM toolbox (http://dreem.openfluor.org/)
% Please refer to the following manuscript and corresponding appendices for more details:

% Murphy K.R., Stedmon C.A., Graeber D. and R. Bro, Fluorescence
%     spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 
%     5, 6557-6566, 2013. DOI:10.1039/c3ay41160e. 

%% Configuration

clc         % Clears the command window
clear       % Clears the workspace
close all   % Closes any other MATLAB windows (e.g., Figure pop-ups)

% List the location of your data, sample log, and where you want the processed data to be saved
%Location_RawData='C:\Users\Administrator\Documents\MATLAB\TEnvR\TestData_EEM\TestData_Aqualog'; 
%Location_SampleLog='C:\Users\Administrator\Documents\MATLAB\TEnvR\TestData_EEM\SampleLog_Aqualog.csv';

Location_RawData='C:\Users\Administrator\Documents\MATLAB\TEnvR 2021\TestData_EEM\TestData_Aqualog'; 
Location_SampleLog='C:\Users\Administrator\Documents\MATLAB\TEnvR 2021\TestData_EEM\SampleLog_Aqualog.csv';

% The exported files will be saved in this folder:
%Location_ProcessedData='C:\Users\Administrator\Documents\MATLAB\TEnvR\TestData_EEM\TestData_Aqualog_Processed'; 

Location_ProcessedData='C:\Users\Administrator\Documents\MATLAB\TEnvR 2021\TestData_EEM\TestData_Aqualog_Processed';
%% Stage 1: Import EEM data files, creates a cube "X" with all the data

cd(Location_RawData)
[X,Emmat,Exmat,filelist_eem,~]=readineems(3,'dat','A4..IC128',[0 1],0,0);
Ex=Exmat(1,:); % Since all files have the same excitation wavelengths
Em=Emmat(:,1); % Since all files have the same emission wavelengths

%% Stage 2: Align the imported data with the sample log

cd(Location_ProcessedData) % Changes the directory

% Start a diary - this function will record everything from the Command
% Window and export a file in the end for keeping a record of the processing
diary EEM_DataCorrectionLog.txt

% Import the sample log
% In the square braket, specify which columns contain text (1) and which numbers (0)
SampleLog=readlogfile(Location_SampleLog,[0 1 1 1 0 1 1 1 1 1 1 1 1 1 0 1]);

% Align the imported EEM data files with according to the SampleLog. 
AnalDate=   alignds(SampleLog,{'EEMfile',filelist_eem},{'AnalDate'});
Cruise=     alignds(SampleLog,{'EEMfile',filelist_eem},{'Cruise'}); 
Site=       alignds(SampleLog,{'EEMfile',filelist_eem},{'Site'}); 
Replicates= alignds(SampleLog,{'EEMfile',filelist_eem},{'Rep'});
SampleID=   alignds(SampleLog,{'EEMfile',filelist_eem},{'SampID'});
QS_Slope=   alignds(SampleLog,{'EEMfile',filelist_eem},{'QS_Slope'});
DF=         alignds(SampleLog,{'EEMfile',filelist_eem},{'dilutionfactor'});

% Turn X into a dataset structure, DS;
DS=assembledataset(X,Ex,Em,'AU','filelist',filelist_eem,[]);

%% Stage 3: Correct the EEMs

% As most of the correction steps have been already performed by the instrument's software prior data export, 
% there is no need to use FDOMcorrect for corrections. Data only need to be undiluted.

% Undilution - rescaling the EEMs corresponding to the dilution factor (DF)
DS_undiluted=DS;
Datacube_DS_undiluted=undilute(DS.X,DF);
DS_undiluted.X=Datacube_DS_undiluted;

%% Stage 4: EEM evaluation and fine correction

% View the EEMs and evaluate if any further corrections are needed
eemview(DS_undiluted,'X',[],[],[],[],[],'rotate','colorbar',[],[])

% Remove residual scattering
% Ray1 = the value used for Rayleigh masking by the instrument's software (10 nm)
% Ray2 = 2 x Ray1 = 20 nm
% Ram1 and Ram2 = 10 nm (the value used for Rayleigh masking by the instrument's software)
DS_Smooth=smootheem(DS_undiluted,[10 10],[10 10],[20 20],[10 10],[0 0 0 0],[],3382,'pause');

% EEMs can be further corrected by removing noise (wavelenghts at which there is no signifficant signal)
DS_Denoised=subdataset(DS_Smooth,[],DS_Smooth.Em>650,DS_Smooth.Ex>600);
DS_Denoised=subdataset(DS_Denoised,[],DS_Denoised.Em<300,DS_Denoised.Ex>500);

% View the EEMs again for quality control after the fine corrections above
eemview(DS_Denoised,'X',[],[],[],[],[],'rotate','colorbar',[],[])

% It was observed that there is an instrument error on EEM 11 - remove using zap
DS_Denoised=zap(DS_Denoised,11,[535 555],[238 242]);

% View the EEMs one last time to verify they do not need any further refinement
eemview(DS_Denoised,'X',[],[],[],[],[],'rotate','colorbar',[],[])

%% Stage 5: Export processed data as csv files

Xout=DS_Denoised; % X out is the dataset to be exported. 

for i=1:Xout.nSample
filename=deblank(char(Xout.filelist{i}));
filename=filename(1:end-4);                 % This removes the extension from filename(.dat)
eem_i=squeeze(Xout.X(i,:,:));               % This converts the 3D dataset to an exportable 2D matrix
eem_i=[[NaN; Xout.Em] [Xout.Ex'; eem_i]];   % This attaches excitation and emission wavelengths to the EEM matrix
csvwrite([Location_ProcessedData '\' filename '_Final.csv'],eem_i) % Export as a .csv file
end

save EEM_DataProcessing_Aqualog.mat         % Exports all variables from the Workspace in a .mat file
diary off

disp('EEM_Process_Agualog - Complete!')