function UVVIS_Process(filename_sample,filename_blank,DF)
%% Description

% This code takes a UV-VIS spectrum (wavelength in column one, absorbance values in column two)
% and converts raw instrumental output into processed data. 
% Blank spectrum can be also loaded (optional). 
% The spectra must be located in .csv files (comma-separated values). 

% Examples:
% UVVIS_Process('Sample 1.csv','Blank.csv',0.0909)
% UVVIS_Process('Sample 1.csv','NoBlank',0.0909)

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

%% Configuration

close all
clearvars -except filename_sample filename_blank DF  

% Parameters that can be modified
DF_blank=1;         % Dilution factor for the blank 
Pathlength=1;       % Cuvette pathlength, in cm
Correction='yes';   % To enable scattering correction = 'yes', to disable = 'no'

% If you want to disable the blank-correction, use 'NoBlank' as filename_blank

%% Stage 1: Load data and interpolate if necessary

Data_Sample=readmatrix(filename_sample); 

if strcmp(filename_blank,'NoBlank') == 1
    blank_correction = 0;
else
    blank_correction = 1;
end

if blank_correction
    Data_Blank=readmatrix(filename_blank);
else
    Data_Blank=[Data_Sample(:,1),zeros(size(Data_Sample(:,1),1),1)];
end

% Sort
Data_Sample_Sorted = sortrows(Data_Sample,1); % Sorts data short -> long wavelength
Data_Blank_Sorted = sortrows(Data_Blank,1);   % Sorts data short -> long wavelength

% Extract data
Wavelengths_Sample=Data_Sample_Sorted(:,1);
Wavelenghts_Blank=Data_Blank_Sorted(:,1);

Absorbance_Sample=Data_Sample_Sorted(:,2);
Absorbance_Blank=Data_Blank_Sorted(:,2);

% Interpolation if wavelengths are not spaced by 1 nm
if Wavelengths_Sample(2)-Wavelengths_Sample(1) ~=1
    
    Wavelengths_Sample_Interpolated=transpose(round(Wavelengths_Sample(1)):1:round(Wavelengths_Sample(end)));
    Wavelengths_Blank_Interpolated=transpose(round(Wavelenghts_Blank(1)):1:round(Wavelenghts_Blank(end)));
    
    Absorbance_Sample_Interpolated=interp1(Wavelengths_Sample,Absorbance_Sample,Wavelengths_Sample_Interpolated);
    Absorbance_Blank_Interpolated=interp1(Wavelenghts_Blank,Absorbance_Blank,Wavelengths_Blank_Interpolated);
    
    Absorbance_Sample=Absorbance_Sample_Interpolated;
    Absorbance_Blank=Absorbance_Blank_Interpolated;
    Wavelengths_Sample=Wavelengths_Sample_Interpolated;
    Wavelenghts_Blank=Wavelengths_Blank_Interpolated;
end

% Interpolation of blank if not of same size as the sample
if size(Wavelengths_Sample,1) ~= size(Wavelenghts_Blank,1)

    Absorbance_Blank_Interpolated=interp1(Wavelenghts_Blank,Absorbance_Blank,Wavelengths_Sample);

    Absorbance_Blank=Absorbance_Blank_Interpolated;
    Wavelenghts_Blank=Wavelengths_Sample;
end

%% Stage 2: Corrections

Data_Sample_New=[Wavelengths_Sample,Absorbance_Sample];
Data_Blank_New=[Wavelenghts_Blank,Absorbance_Blank];

% Undilute the sample
Data_Sample_Undiluted=[Data_Sample_New(:,1),Data_Sample_New(:,2)/DF];
Data_Blank_Undiluted=[Data_Blank_New(:,1),Data_Blank_New(:,2)/DF_blank];

% Blank-correct the sample
Data_Sample_Undiluted_Blanked=[Data_Sample_Undiluted(:,1),(Data_Sample_Undiluted(:,2)-Data_Blank_Undiluted(:,2))];

Wavelength_end=Data_Sample_Undiluted_Blanked(end,1);
if Wavelength_end <= 0
    Correction = 'no';
elseif Wavelength_end >= 800
    Wavelength_end = 800;
end

% Scattering correction
if strcmp(Correction,'yes')
    if Wavelength_end >= 701
        [Abs700800,~]=Absorbance_Range(Data_Sample_Undiluted_Blanked,701,Wavelength_end);
        Abs700800_AVG=mean(Abs700800);
        Data_Sample_Undiluted_Blanked_Smooth=[Data_Sample_Undiluted_Blanked(:,1),(Data_Sample_Undiluted_Blanked(:,2)-Abs700800_AVG)];
    else
        Data_Sample_Undiluted_Blanked_Smooth=Data_Sample_Undiluted_Blanked;
    end
else
    Data_Sample_Undiluted_Blanked_Smooth=Data_Sample_Undiluted_Blanked;
end

% Correct for cuvette pathlength, units now are cm-1
Data_Sample_Undiluted_Blanked_Denoised_Pathcorr=[Data_Sample_Undiluted_Blanked_Smooth(:,1),(Data_Sample_Undiluted_Blanked_Smooth(:,2)/Pathlength)];

% Final data
Data_Final_Decadic=Data_Sample_Undiluted_Blanked_Denoised_Pathcorr;

%% Stage 3: Export

dlmwrite([filename_sample(1:end-4) '_Final.csv'],Data_Final_Decadic)
disp(['Finished applying UVVIS_Process on ' char(filename_sample) ' (blank:' char(filename_blank) ')'])
end

%% Internal Functions

% Extracts absorbance values over a range of wavelengths 
function [Abs_Range,Wav_Range] = Absorbance_Range(Data,Start,End)

Wav_Range=Start:1:End;

for j=1:size(Wav_Range,2)
    Abs_Range(j)=Absorbance_OnePoint(Data,Wav_Range(j));
end

end

% Extracts an absorbance value at a specified wavelength
function Abs_OnePoint = Absorbance_OnePoint(Data,Wavelength)

[~,row] = ismember(Wavelength,Data(:,1));
Abs_OnePoint = Data(row,2);
end