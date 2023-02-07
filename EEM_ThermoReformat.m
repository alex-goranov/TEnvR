%% Description
% This code takes one or multiple processed EEM spectra output by a 
% Thermo Scientific Lumina spectrfluorometer and converts them into a regular csv format.
% The EEM spectrum must be located in a .xls files (Microsoft Excel 97-2003 Worsheet). 

% Example: 
% Run the code by typing EEM_ThermoReformat in the command windwow.

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

%% Code

close all
clearvars 

Files=dir('*_3DScan.xls');

for z=1:size(Files,1)
    FileName=Files(z).name;
    ThermoConvert(FileName)
end
    
disp('EEM_ThermoConvert - Complete!')

%% Internal Functions

% This function performs the actual conversion 
function ThermoConvert(filename)
dataM=xlsread(filename); % Loads the data from the xls export of the Thermo Fluorometer

Signal=dataM(:,3:3:end); % data is in the 3rd column, and then every consecutve 3rd

Excitation=dataM(1,1:3:end); % data is in the 3rd row, 1st column, and then every consecutive 3rd

Emission=dataM(:,2); % emission wavelengths are in the second column

xlswrite([filename(1:end-11) '_ref.xlsx'], Excitation, 1, 'B1');
xlswrite([filename(1:end-11) '_ref.xlsx'], Emission, 1, 'A2');
xlswrite([filename(1:end-11) '_ref.xlsx'], Signal, 1, 'B2');

disp(['Finished converting ' char(filename)])
end