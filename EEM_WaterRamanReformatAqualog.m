%% Description

% This code takes multiple Water Raman spectra exported by an Aqualog fluorometer
% Filenames must end with " - RSUNIT (01) - Raman Area Graph.dat". 
% Data is exported as two-dimensional arrays in .csv files.

% Example: 
% Run the code by typing EEM_WaterRamanReformatAqualog in the command windwow.

%% Copyright and License Notice: 

% Copyright © 2022 Old Dominion University Research Foundation, Norfolk VA, USA
% All rights reserved.

% This file is part of the Toolbox for Environmental Research (TEnvR). Please cite the toolbox as follows: 
% Goranov, A. I., Sleighter, R. L., Yordanov, D. A., and Hatcher, P. G. (2023): 
% TEnvR: MATLAB-Based Toolbox for Environmental Research, Analytical Methods, doi: 10.1039/d3ay00750b.

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

% Clear up your workspace and close any other windows. 
clear
close all

% Import filenames
Files=dir('*Raman Area Graph.dat');

for i=1:size(Files,1)
    Filename=Files(i).name;
    WaterRamanReformat(Filename);
end

%% Internal Functions

% This function performs the actual reformatting
function WaterRamanReformat(Filename)

data=importdata(Filename);
wavelengths=data.data(:,3); 
intensity=data.data(:,4);
indexes=find(isnan(wavelengths) == 1);
wavelengths(indexes)=[];
intensity(indexes)=[];

dlmwrite([Filename(1:end-23) '_Raman_ref.csv'],[wavelengths,intensity], 'delimiter', ',')

disp(['Finished reformatting ' char(Filename)])
end