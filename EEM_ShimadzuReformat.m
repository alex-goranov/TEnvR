%% Description
% This code takes one or multiple processed EEM spectra output by a 
% Shimadzu RF-6000 spectrfluorometer (.txt files with lots of text)and 
% converts them into a readable format by the drEEM toolbox.

% Example: 
% Run the code by typing EEM_ShimadzuReformat in the command windwow.

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

close all
clearvars 

Files=dir('*.txt'); 

for z=1:size(Files,1)
  FileName=Files(z).name;
  EEM_ShimadzuTXTtoCSV(FileName);
end

disp('EEM_ShimadzuConvert - Complete!')

%% Internal Functions

% This function performs the actual conversion 
function EEM_ShimadzuTXTtoCSV(filename)

Input=readmatrix(filename);
Input(1,1)=0;
writematrix(Input,[filename(1:end-4) '_ref.csv'])

disp(['Finished converting ' char(filename)])
end