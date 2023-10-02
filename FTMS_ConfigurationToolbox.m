classdef FTMS_ConfigurationToolbox
%% Description 

% This is a configuration code for the codes in TEnvR for assessing
% FT-ICR-MS data (comparisons, statistics, figures, etc.). 

% The parameters in this file are not used for data processing by the 
% FTMS_RefinementPeaks, FTMS_FormulaAssignment, FTMS_RefinementFormulas 
% codes - those parameters are found in the FTMS_ConfigurationAssignment code.

% Please record any changes that you make in your laboratory notebook. We
% also recommend to save a separate configuration file for each different
% project. Please note that if you change the filename of the file in the
% TEnvR folder, you also must change the filename on line 1 above
% corresponding to classdef. 

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
    
%% Configuration parameters
    properties
    
    % Define mass accuracy and precision 
    MassAccuracy = 1;   % 1ppm for 12T ODU and 15T OSU FT-ICR-MS instruments
    Precision = 5;      % Precise m/z measurements up to 5 decimals for 12T ODU and 15T OSU FT-ICR-MS instruments. 6th decimal is the uncertainty.
    
    % Define ionization mode and heteroelements 
    Ionization_Mode = 'negative';   % Choose ionization mode: 'negative' or 'positive'
    Heteroelement_Halogen = false;  % Select if the used heteroelement is a halogen or halogen-like (e.g., F, Cl, Br, I, D)
    Heteroelement_N15 = false;      % Select if the used heteroelement is N15 or N-like.
    Mass_Heteroelement = 0;         % Choose mass of heteroelement. Use 0 if none. 35Cl =  34.96885271

    % Define region of spectrum of highest reliability - used for alignment & statitics 
    MZ_low  = 300;                  % Default = 300
    MZ_high = 800;                  % Default = 800 
    
    % Statistical parameters
    CL_alpha = 95;                  % Default = 95 % (equivalent to p-value threshold of 0.05)
    MinSamples = 1;                 % Default = 1
    PresenceAbsence = false;        % Default = false
    pvalue_adjustment = false;       % Default = false
    p_adjust_method = 'BH';         % Default = 'BH'
               
    % File format - describes the format of the _Final.xlsx sheets.
    % DO NOT CHANGE unless 1) altering the format of the _Final.xlsx sheets in FTMS_FormulaRefinement
    %                      2) Planning on using TEnvR on formula lists of different format
    
    Column_mz=1;            Column_C=6;     Column_P=11;            Column_HC=16;
    Column_Magnitude=2;     Column_H=7;     Column_E=12;            Column_DBE=17;
    Column_SN=3;            Column_O=8;     Column_ExactMass=13;    Column_DBEC=18;    
    Column_EstimC=4;        Column_N=9;     Column_error=14;        Column_AImod=19;
    Column_Type=5;          Column_S=10;    Column_OC=15;           Column_RelMagnitude=20; % For FTMS_Metrics and FTMS_Compare
   
    % Elemental masses
    Mass_12C = 12.0000000;      Mass_1H = 1.0078250;        Mass_16O = 15.9949146;
    Mass_14N = 14.0030740;      Mass_32S = 31.9720707;      Mass_31P = 30.9737615;
    Mass_40K = 39.963998700;    Mass_23Na = 22.9897697;     Mass_electron = 0.0005485799; 
    end
end