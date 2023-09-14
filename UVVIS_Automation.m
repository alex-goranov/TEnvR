%% Description

% This code takes multiple UV-VIS spectra (wavelength in column one, absorbance values in column two)
% and automatically processes them, derives them, calculates the metrics for all of them,
% and organizes all of the data in a large Excel file.
% The spectra must be located in .csv files (comma-separated values). 
 
% Example: 
% Run the code by typing UVVIS_Automation in the command windwow.

%% Copyright and License Notice: 

% Copyright © 2022 Old Dominion University Research Foundation, Norfolk VA, USA
% All rights reserved.

% This file is part of the Toolbox for Environmental Research (TEnvR). Please cite the toolbox as follows: 
% Goranov, A. I., Sleighter, R. L., Yordanov, D. A., and Hatcher, P. (2023): 
% TEnvR: MATLAB-Based Toolbox for Environmental Research, Analytical Methods, doi: XXXXXXXXXXX.

% TEnvR is free software for non-commercial use: you can redistribute it and/or modify 
% %it under the terms of the GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version. 
% Users wishing to use TEnvR commercially must obtain a commercial license from Old Dominion University Research Foundation. 

% TEnvR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

% You should have received a copy of the GNU General Public License along with TEnvR
% (located in the Supplementary files directory of the toolbox). If not, see <https://www.gnu.org/licenses/>.

% Please contact the corresponding authors Drs. Aleksandar Goranov (aleksandar.i.goranov@gmail.com) 
% and Patrick Hatcher (phatcher@odu.edu) with any questions or concerns.

%% Data Import & Configuration

% Clear up your workspace and close any other windows. 
clear
close all

% You will need to key in the name of the file containing the UV-VIS
% spectrum of the blank. Input without apostrophes. E.g.: Blank.csv
blankname=inputdlg('Enter the filename of the blank','Blank',1);

% Import filenames
Files=dir('*.csv');
Files=natsortfiles(Files);

% Remove the blank from the structure containing the spectra of the actual samples
Files_Filenames={Files.name};
index = find(strcmp(Files_Filenames,blankname));
Files(index)=[];

answer = questdlg('Were any of the samples diluted?', ...
	'Dilution', ...
	'Yes','No','No');
% Handle response
switch answer
    case 'Yes' %Import dilution factors that need to be manually keyed in
        for i=1:size(Files,1)
            FileName=Files(i).name;
            DF_text(i,:)=inputdlg(['Enter dilution factor for sample "' FileName '"'],'DF import',1);
        end
        DF=str2double(DF_text);
    case 'No' %All DF = 1
    DF=zeros(1,size(Files,1))+1;
    DF=DF';
end

%% Apply UVVIS_Process.m on all samples

for z=1:size(Files,1)
    FileName=Files(z).name;
    UVVIS_Process(FileName,string(blankname),DF(z))
end

%% Apply UVVIS_Derivative.m on all samples

Files_final=dir('*_Final.csv');
Files_final=natsortfiles(Files_final);
for i=1:size(Files_final,1)
    FileName=Files_final(i).name;
    D=UVVIS_Derivative(FileName);
    Matrix_wavelspectra(:,i)=D(:,1);
    Matrix_derivspectra(:,i)=D(:,2);
end

%% Apply UVVIS_Metrics.m on all samples

for i=1:size(Files_final,1)
    FileName=Files_final(i).name;
    [M1,M2,M3]=UVVIS_Metrics(FileName);
    Matrix_metrics(:,i)=M1;
    Matrix_decspectra(:,i)=M2(:,2);
    Matrix_Napspectra(:,i)=M3(:,2);
end

Matrix_norm=Matrix_decspectra;
Wavelengths_norm=M2(:,1);

% replace negative and NaN values with zero
Matrix_norm(Matrix_norm<0)=0; 
Matrix_norm(isnan(Matrix_norm))=0;

% Remove variables above 450 nm.
indexes=find(Wavelengths_norm > 450);

Matrix_norm(indexes,:)=[];
Wavelengths_norm(indexes)=[];

% Normalize 
TotalSpectralInt=sum(Matrix_norm,1);
for i=1:size(Matrix_norm,2)
    RawInt=Matrix_norm(:,i);
    Norm(:,i)=RawInt./TotalSpectralInt(i);
end

%% Export a Master Report

warning('off','MATLAB:xlswrite:AddSheet')
labels={'Abs-230 (dec)' 'Abs-254 (dec)' 'Abs-280 (dec)' 'Alpha350' 'Alpha325' 'Alpha300' 'Alpha250' 'Alpha254' 'E4toE6' 'E2toE3' 'Total CDOM' 'Slope(275-295)' 'Slope(290-320)' 'Slope(350-400)' 'Slope(350-450)' 'SR'};
Files_final_Filenames={Files_final.name};
Files_final_Filenames = cellfun(@(x) x(1:end-4), Files_final_Filenames, 'un', 0);

% Export Decadic spectra in Sheet1
xlswrite('UVVIS_MasterReport.xlsx',Files_final_Filenames,1,'B1')
xlswrite('UVVIS_MasterReport.xlsx',Matrix_decspectra,1,'B2')
xlswrite('UVVIS_MasterReport.xlsx',M2(:,1),1,'A2')

xlswrite('UVVIS_MasterReport.xlsx',Files_final_Filenames,'Normalized','B1')
xlswrite('UVVIS_MasterReport.xlsx',Norm,'Normalized','B2')
xlswrite('UVVIS_MasterReport.xlsx',Wavelengths_norm,'Normalized','A2')

% Export Napierian spectra in sheet "Napierian"
xlswrite('UVVIS_MasterReport.xlsx',Files_final_Filenames,'Napierian','B1')
xlswrite('UVVIS_MasterReport.xlsx',Matrix_Napspectra,'Napierian','B2')
xlswrite('UVVIS_MasterReport.xlsx',M3(:,1),'Napierian','A2')

% Export First derivative in sheet "1st Deriv"
xlswrite('UVVIS_MasterReport.xlsx',Files_final_Filenames,'1st Deriv','B1')
xlswrite('UVVIS_MasterReport.xlsx',Matrix_derivspectra,'1st Deriv','B2')
xlswrite('UVVIS_MasterReport.xlsx',Matrix_wavelspectra(:,1),'1st Deriv','A2')

% Export all Metrics in sheet "Metrics"
xlswrite('UVVIS_MasterReport.xlsx',labels,'Metrics','C1')
xlswrite('UVVIS_MasterReport.xlsx',Matrix_metrics','Metrics','C2')
xlswrite('UVVIS_MasterReport.xlsx',{'DF'},'Metrics','B1')
xlswrite('UVVIS_MasterReport.xlsx',DF,'Metrics','B2')
xlswrite('UVVIS_MasterReport.xlsx',Files_final_Filenames','Metrics','A2')

disp('UVVIS_Automation - Complete!')