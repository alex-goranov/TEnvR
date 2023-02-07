%% Description

% This code takes multiple UV-VIS spectra (wavelength in column one, absorbance values in column two)
% and calculates differential ("difference") spectra. The files must have
% already been processed and thus labeled as "_Final". 
% All spectra must be located in .csv files (comma-separated values). 

% Example: 
% Run the code by typing UVVIS_Differentiation in the command windwow.

%% Copyright and License Notice: 

% Copyright © 2022 Old Dominion University Research Foundation, Norfolk VA, USA
% All rights reserved.

% This file is part of the Toolbox for Environmental Research (TEnvR). Please cite the toolbox as follows: 
% Goranov, A. I., Sleighter, R. L., Yordanov, D. A., and Hatcher, P. (2023): 
% TEnvR: MATLAB-Based Toolbox for Environmental Research, Journal TBD, doi: XXXXXXXXXXX.

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

%% Stage 1: Data Import & Configuration

% Clear up your workspace and close any other windows. 
clear
close all

% You will need to key in the name of the sample which will be subtracted
% from all other files (e.g., a control sample)
% Input without apostrophes. E.g.: Sample 1_Final.csv
controlname=inputdlg('Enter the filename of the control','Control',1);

% Import filenames
Files=dir('*_Final.csv');
Files=natsortfiles(Files);

% Remove the control from the structure containing the rest of the processed spectra
Files_Filenames={Files.name};
index = find(strcmp(Files_Filenames,controlname));
Files(index)=[];

if index == 0
    error('Error! The sample name you put in for the control does not exist!')
end

%% Stage 2: Differentiate all spectra
Data_control=readmatrix(char(controlname)); 
Data_control=sortrows(Data_control,1);

for i=1:size(Files,1)
    FileName=Files(i).name;
    Data_sample=readmatrix(FileName); 
    Data_sample=sortrows(Data_sample,1);
    
    if size(Data_sample(:,1),1) ~= size(Data_control(:,1),1)
    error('Error! Sample and Blank have different wavelength ranges!')
    end
    
    DifferenceSpectrum=Data_sample(:,2)-Data_control(:,2);
    Matrix_derivspectra(:,i)=DifferenceSpectrum;
    Matrix_originalspectra(:,i)=Data_sample(:,2);
    disp(['Finished differentiating ' char(FileName)])
end

%% Stage 3: Export a Master Report
warning('off','MATLAB:xlswrite:AddSheet')

% Exporting original spectra to Sheet1
xlswrite('UVVIS_Differentiation.xlsx',{Files.name},1,'C1')
xlswrite('UVVIS_Differentiation.xlsx',controlname,1,'B1')
xlswrite('UVVIS_Differentiation.xlsx',Data_control,1,'A2')
xlswrite('UVVIS_Differentiation.xlsx',Matrix_originalspectra,1,'C2')

% Exporting differential spectra to sheet "Differential"
xlswrite('UVVIS_Differentiation.xlsx',{Files.name},'Differential','C1')
xlswrite('UVVIS_Differentiation.xlsx',strcat(controlname,' CONTROL'),'Differential','B1')
xlswrite('UVVIS_Differentiation.xlsx',Data_control,'Differential','A2')
xlswrite('UVVIS_Differentiation.xlsx',Matrix_derivspectra,'Differential','C2')

disp('UVVIS_Differentiation - Complete!')