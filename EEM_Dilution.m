function EEM_Dilution (filename,action,Vin,Vfinal)
%% Description

% This code takes an EEM spectrum and scales up (undilutes) or scales down (dilutes)
% the intensity data based on a dilution factor (Vin=volume you took and diluted to Vfin). 
% The spectrum must be located in a .csv file (comma-separated values). 

% Can use this code to normalize your absorbance to an external parameter 
% such as dissolved organic carbon DOC content:

% Examples:
% EEM_Dilution('Sample 1_ref_Final.csv','dilute',1,5)
% EEM_Dilution('Sample 1_ref_Final.csv','undilute',1,10)
% EEM_Dilution('Sample 1_ref_Final.csv','undilute',DOC,1)
% EEM_Dilution('Sample 1_ref_Final_undil_2.csv','dilute',DOC,1)

% In the examples above: (Vin = 1 mL; Vfin = 5 mL or 10 mL, DOC = 2 ppm-C)

%% Copyright and Cicense Notice: 

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

Data=csvread(filename); % Import data
Ex=Data(1,:); Em=Data(:,1); % Extract excitation and emission ranges
Ex=Ex(2:end); Em=Em(2:end); % Extract excitation and emission ranges
Int=Data; Int(1,:)=[]; Int(:,1)=[]; % Extract intensity data 
DF=Vin/Vfinal;

% dilute or undilute
if action == "dilute"
    Int_dil=Int*DF;
    filename_new=[char(filename(1:end-4)) '_dil_' num2str(round(DF,4)) '.csv'];
    dlmwrite(filename_new,[NaN,Ex;[Em,Int_dil]], 'delimiter', ',', 'precision', '%i');
elseif action == "undilute"
    Int_dil=Int/DF;
    filename_new=[char(filename(1:end-4)) '_undil_' num2str(round(DF,4)) '.csv'];
    dlmwrite(filename_new,[NaN,Ex;[Em,Int_dil]], 'delimiter', ',', 'precision', '%i')
    
else
    error('Error! Your "action" must be either "dilute" or undilute"')
end
disp(['Finished applying EEM_dilution on ' char(filename)])
end