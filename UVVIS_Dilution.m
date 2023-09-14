function UVVIS_Dilution (filename,action,Vin,Vfinal)
%% Description

% This code takes a UV-VIS spectrum (wavelength in column one, absorbance values in column two)
% and scales up (undilutes) or scales down (dilutes) the absorbance data based on a dilution factor
% (Vin=volume you took and diluted to Vfin). The spectrum must be located
% in a .csv file (comma-separated values). 

% Can use this code to normalize your absorbance to an external parameter 
% such as dissolved organic carbon (DOC) content.

% Examples:
% UVVIS_Dilution('UVVIS_TestSample1.csv','dilute',4,16)
% UVVIS_Dilution('UVVIS_TestSample1_dil_0.25.csv','undilute',4,16)
% UVVIS_dilution ('filename','undilute',23,1) 

% In the examples above: (Vin = 4 mL; Vfin = 16 mL, DOC = 23 ppm-C)

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

%% Data Import 

data=dlmread(filename);
data=sortrows(data,1); % Sorts data short -> long wavelength
wavelengths=data(:,1);
absorbance=data(:,2);
DF=Vin/Vfinal;

filename_orig=filename;

if action == "dilute"
    absorbance_dil=absorbance*DF;
    filename=[char(filename(1:end-4)) '_dil_' num2str(round(DF,4)) '.csv'];
    dlmwrite(filename,[wavelengths,absorbance_dil], 'delimiter', ',', 'precision', '%i');
elseif action == "undilute"
    absorbance_undil=absorbance/DF;
    filename=[char(filename(1:end-4)) '_undil_' num2str(round(DF,4)) '.csv'];
    dlmwrite(filename,[wavelengths,absorbance_undil], 'delimiter', ',', 'precision', '%i')
    
else
    error('Error! Your "action" must be either "dilute" or undilute"')
end

disp(['Finished applying UVVIS_dilution on ' char(filename_orig)])
end