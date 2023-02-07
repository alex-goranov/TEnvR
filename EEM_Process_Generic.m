%% Description

% This code is used to process raw EEM spectra that have been acquired on  
% spectrofluorometers other than a Horiba Aqualog. 
% We show an example here with data from Shimadzu RF-6000.

% The spectra have been only corrected for Instrument-specific responces 
% (automatically done by the software)

% Spectra for this version of the code must have been exported as .csv files.
% Note that the "readineems" code of drEEM is capable of importing EEM 
% files in various formats: 'csv', 'xls', 'xlsx','dat', 'txt' - modify
% below if needed. 

% The EEM_Process_Generic script is based on the tutorials by Murphy et al. 2013. 
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
Location_RawData='C:\Users\Administrator\Documents\MATLAB\TEnvR 2021\TestData_EEM\TestData_Generic'; 
Location_SampleLog='C:\Users\Administrator\Documents\MATLAB\TEnvR 2021\TestData_EEM\SampleLog_Generic.csv';

% The exported files will be saved in this folder
Location_ProcessedData='C:\Users\Administrator\Documents\MATLAB\TEnvR 2021\TestData_EEM\TestData_Generic_Processed';

%% Stage 1: Import EEM data and other auxiliary files, creates a cube "X" with all the EEM data

cd([Location_RawData '\BlankEEMs']) 
[X_b,Emmat_b,Exmat_b,filelist_b,outdata_b]=readineems(1,'csv','A1..BD96',[1 1],0,0);
Exb=Exmat_b(1,:); % Since all files have the same excitation wavelengths
Emb=Emmat_b(:,1); % Since all files have the same emission wavelengths

cd([Location_RawData '\CorrectionFiles'])
Excor=csvread('ExcitationCorr.csv');
Emcor=csvread('EmissionCorr.csv');

cd([Location_RawData '\RawEEMs'])
[X,Emmat,Exmat,filelist_eem,outdata]=readineems(1,'csv','A1..BD96',[1 1],0,0); 
Ex=Exmat(1,:); % Since all files have the same excitation wavelengths
Em=Emmat(:,1); % Since all files have the same emission wavelengths

cd([Location_RawData '\UVVIS_IFE'])
[S_abs,W_abs,wave_abs,filelist_abs]=readinscans('Abs','csv','A1..B571',0,0);

cd([Location_RawData '\WaterRaman350'])
[S_R,W_R,wave_R,filelist_R]=readinscans('R350','csv','A1..B171',0,0);
RamEx=350;     % Excitation wavelength of water

%% Stage 2: Align the imported data with the sample log

cd(Location_ProcessedData) % Changes the directory

% Start a diary - this function will record everything from the Command
% Window and export a file in the end for keeping a record of the processing
diary EEM_DataCorrectionLog.txt

% Import the sample log
% In the square braket, specify which columns contain text (1) and which numbers (0)
SampleLog=readlogfile(Location_SampleLog,[0 1 1 1 0 1 1 1 1 1 1 1 1 1 0 1]);

% Align the imported EEM data files with according to the SampleLog. 
dates=      alignds(SampleLog,{'EEMfile',filelist_eem},{'AnalDate'});
cruises=    alignds(SampleLog,{'EEMfile',filelist_eem},{'Cruise'}); 
sites=      alignds(SampleLog,{'EEMfile',filelist_eem},{'Site'}); 
replicates= alignds(SampleLog,{'EEMfile',filelist_eem},{'Rep'});
sampleID=   alignds(SampleLog,{'EEMfile',filelist_eem},{'SampID'});
Q=          alignds(SampleLog,{'EEMfile',filelist_eem},{'QS_Slope'});
DF=     alignds(SampleLog,{'EEMfile',filelist_eem},{'dilutionfactor'});

% Align the imported EEM data files with the other imported data (blank files, UV-VIS files, Water Raman files)
Sabs=       alignds(SampleLog,{'EEMfile',filelist_eem},{'ABSfile',filelist_abs,S_abs}); 
B=          alignds(SampleLog,{'EEMfile',filelist_eem},{'BlankFile',filelist_b,X_b}); 
Sr=         alignds(SampleLog,{'EEMfile',filelist_eem},{'RamanFile',filelist_R,S_R});

% Turn X into a dataset structure, DS;
DS=assembledataset(X,Ex,Em,'RU','site',sites,'rep',replicates,'filelist',filelist_eem,'ID',sampleID,'cruise',cruises,'date',dates,[]); 

%% Stage 3: Correct the EEMs

% Prepare a Water Raman data matrix (RamMat) of the imported Water Raman spectra
% by attaching emission wavelenegths header to the matrix with the water raman intensity values
RamMat=[wave_R;S_R];  % 2D Matrix of unmatched Raman scans (corresponding to filelist_R)

% Determine the Raman Integration Range following Murphy 2011. 
[IR,IRmed,IRdiff] = ramanintegrationrange(RamMat,filelist_R,RamEx,1800,6,0,0.1);

% Create a variable RamOpt with the obtained parameters for Water Raman normalization
RamOpt=[RamEx IR(1,:)]; 

% Prepare an absorbance data matrix (A) by attaching absorbance wavelenegths
% header to the matrix with the absorbance values
A=[wave_abs;S_abs]; 

% Prepare a Water Raman data matrix (RamMat) of Water Raman spectra (aligned
% with the EEM sample names they correspond to) by attaching emission wavelenegths 
% header to the matrix with the water raman intensity values
W=[wave_R;Sr]; % W = 2D Matrix of matched (350 nm) Raman scans

% Correct the EEM spectra following Murphy et al., 2010  
[XcRU, Arp, IFCmat, BcRU, XcQS, QS_RU]=fdomcorrect(DS,DS.Ex,DS.Em,Emcor,Excor,W,RamOpt,A,B,[],[],[]);

% Undilution - rescaling the EEMs corresponding to the dilution factor (DF)
DS_undiluted=DS;
Datacube_DS_undiluted=undilute(DS.X,DF);
DS_undiluted.X=Datacube_DS_undiluted;

%% Stage 4: EEM evaluation and fine correction

% View the EEMs and evaluate if any further corrections are needed
eemview(DS_undiluted,'X',[],[],[],[],[],'rotate','colorbar',[],[])

% Remove residual scattering
DS_Smooth=smootheem(DS_undiluted,[15 300],[15 15],[350 20],[10 10],[0 0 0 0],[],3382,'pause');

% EEMs can be further corrected by removing noise (wavelenghts at which there is no signifficant signal)
DS_Denoised=subdataset(DS_Smooth,[],DS_Smooth.Em>600,DS_Smooth.Ex>450);
DS_Denoised=subdataset(DS_Denoised,[],DS_Denoised.Em<270,[]);

% View the EEMs again for quality control after the fine corrections above
eemview(DS_Denoised,'X',[],[],[],[],[],'rotate','colorbar',[],[])

%% Stage 5: Export processed data as csv files
Xout=DS_Denoised; % Xout is the dataset to be exported. 

for i=1:Xout.nSample
filename=deblank(char(Xout.filelist{i}));
filename=filename(1:end-4);                 % This removes the extension from filename(.csv)
eem_i=squeeze(Xout.X(i,:,:));               % This converts the 3D dataset to an exportable 2D matrix
eem_i=[[NaN; Xout.Em] [Xout.Ex'; eem_i]];   % This attaches excitation and emission wavelengths to the EEM matrix
csvwrite([Location_ProcessedData '\' filename '_Final.csv'],eem_i) % Export as a .csv file
end

save EEM_DataProcessing_Generic.mat         % Exports all variables from the Workspace in a .mat file
diary off

disp('EEM_Process_Generic - Complete!')