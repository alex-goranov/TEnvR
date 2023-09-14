%% Description

% This code takes one or multiple processed EEM spectra and transposes them.
% The EEM spectrum must be located in a .csv file (comma-separated values). 

% Example: 
% Run the code by typing EEM_Transposition in the command windwow.

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

%% Code

close all
clearvars 

Files=dir('*_Final.csv');
for z=1:size(Files,1)
    FileName=Files(z).name;
    TransposeEEM(FileName)
end

disp('EEM_Transposition - Complete!')

%% Internal Functions

% This function performs the actual transposition 
function TransposeEEM(filename)

Data=csvread(filename);
Ex=Data(1,:); Em=Data(:,1); 
Em=Em(2:end); 
Int=Data;

Int(1,:)=[]; Int(:,1)=[]; 
Int=Int';
Em=Em';
Ex=Ex';

NewMatrix=[Em;Int];
NewMatrix=[Ex,NewMatrix];

csvwrite([filename(1:end-4) '_transposed.csv'],NewMatrix)

disp(['Finished transposing ' char(filename)])
end