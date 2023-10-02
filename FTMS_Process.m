function FTMS_Process(sample,blank,config)
%% Description

% This code takes an FT-ICR-MS file with calibrated data for a sample and a blank and:
    %1) Refines the peak list of the sample using FTMS_RefinementPeaks
    %2) Assigns molecular formulas to the refined peak list using FTMS_FormulaAssignment
    %3) Refined the molecular formulas using FTMS_RefinementFormulas

% Example
% FTMS_Process('Sample 2.txt','Blank PPL.txt',FTMS_ConfigurationAssignment)

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

%% Stage 1: Peak Refinement

FTMS_RefinementPeaks(sample,blank,config)

%% Stage 2: Formula assignment

FTMS_FormulaAssignment([sample(1:end-4) '_Refinement.txt'],config)

%% Stage 3: Formula Refinement

FTMS_RefinementFormulas([sample(1:end-4) '_Refinement_F.txt'],config)

end