%% Description

% This code takes multiple UV-VIS spectra exported by an Aqualog Fluorometer.
% Filename must end with "- Abs Spectra Graphs.dat". The code reformats files 
% into two-dimensional arrays and export them into .csv files.

% Example: 
% Run the code by typing UVVIS_ReformatAgualog in the command windwow.

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

% Clear up your workspace and close any other windows. 
clear
close all

% Import filenames
Files=dir('* Abs Spectra Graphs.dat');

for i=1:size(Files,1)
    Filename=Files(i).name;
    ReformatAqualog(Filename);
end

%% Internal Functions

% Actual function doing the reformatting
function ReformatAqualog(Filename)
data=importdata(Filename);
wavelengths=data.data(:,1); 
absorbance=data.data(:,10);

dlmwrite([Filename(1:end-30) '.csv'],[wavelengths,absorbance], 'delimiter', ',')
disp(['Finished reformatting ' char(Filename)])
end