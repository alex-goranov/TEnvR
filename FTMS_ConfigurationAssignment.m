classdef FTMS_ConfigurationAssignment
%% Description 

% This is a configuration code for the codes performing FT-ICR-MS data
% processing (FTMS_RefinementPeaks, FTMS_FormulaAssignment, FTMS_RefinementFormulas).

% Please record any changes that you make in your laboratory notebook. We
% also recommend to save a separate configuration file for each different
% project. Please note that if you change the filename of the file in the
% TEnvR folder, you also must change the filename on line 1 above
% corresponding to classdef. 

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
    
%% Configuration parameters
    properties
        
    %% PEAK REFINEMENT

    % Define m/z cutoffs
    MZcutoff_low=0;     % Low m/z value
    MZcutoff_high=2000; % High m/z value

    % Export peaks which have associated 37Cl isotopologues (i.e., organochlorine formulas)
    Export_Cl=false;

    % Define mass accuracy and precision 
    MassAccuracy = 1;   % 1ppm for 12T ODU and 15T OSU FT-ICR-MS instruments
    Precision = 5;      % Precise m/z measurements up to 5 decimals for 12T ODU and 15T OSU FT-ICR-MS instruments. 6th decimal is the uncertainty.

    %% FORMULA ASSIGNMENT
    
    Ionization_Mode = 'negative';   % Choose ionization mode: 'negative' or 'positive'
    Heteroelement_Halogen = false;  % Select if the used heteroelement is a halogen or halogen-like (e.g., F, Cl, Br, I, D)
    Heteroelement_N15 = false;      % Select if the used heteroelement is N15 or N-like.
    Mass_Heteroelement = 0;         % Choose mass of heteroelement. Use 0 if none. 35Cl =  34.96885271

    Mass_12C = 12.0000000;      Mass_1H = 1.0078250;        Mass_16O = 15.9949146;
    Mass_14N = 14.0030740;      Mass_32S = 31.9720707;      Mass_31P = 30.9737615;
    Mass_40K = 39.963998700;    Mass_23Na = 22.9897697;     Mass_electron = 0.0005485799;

    % Elemental criteria (ranges) 
    C_min = 5;      H_min = 5;      O_min = 1;      N_min = 0;
                    H_max = 100;    O_max = 30;     N_max = 5;

    S_min = 0;      P_min = 0;      E_min = 0;      K_min = 0;      Na_min = 0;
    S_max = 4;      P_max = 2;      E_max = 0;      K_max = 0;      Na_max = 0;

    %% FORMULA REFINEMENT using Elemental Constrains (mainly based on Stubbins et al., 2010)

    %Stubbins, A., Spencer, R.G.M., Chen, H.M., Hatcher, P.G., Mopper, K., Hernes, P.J., Mwamba, V.L.,
    %Mangangu, A.M., Wabakanghanzi, J.N. and Six, J. (2010) Illuminated darkness: Molecular signatures 
    %of Congo River dissolved organic matter and its photochemical alteration as revealed by ultrahigh 
    %precision mass spectrometry. Limnology and Oceanography 55, 1467-1477.

    O_to_C_min = 0;             H_to_C_min = 0.33;      N_to_C_min = 0;
    O_to_C_max = 1.2;           H_to_C_max = 2.25;      N_to_C_max = 0.5;
    O_to_C_delta_max = 2;

    S_to_C_min = 0;             P_to_C_min = 0;         O_to_P_min = 3;
    S_to_C_max = 0.2;           P_to_C_max = 0.1;       O_to_P_max = Inf;

    DBE_min = 0;
    DBE_max = 50; 

    %% FORMULA REFINEMENT using filters 

    % Stage 3: Isotopic Filter - fitler formulas based on C-number
    Filter_Isotopes = true;             % Enable or disable filter (true or false)
    Filter_Isotopes_Range = 1.0;        % Allowed range for the estimated C-number (+/-)
    Filter_Isotopes_SNthreshold = 25;   % Minimum S/N for the 13C isotopologue peak to be used for C-estimation

    % Stage 3: Uniqueness Filter - fitler formulas based on ambiguity
    Filter_Unique = true; %Enable or disable filter (true or false)

    % Stage 4: KMD filter - fitler formulas based KMD homologous series
    Filter_KMD = true; % Enable or disable filter (true or false)      
    Filter_KMD_Series = double([14/14.015650, 2/2.015650, 44/43.989829, ...
        30/30.010565, 44/44.026215, 16/15.994915, 18/18.010565, 17/17.026549, ...
        36/35.976678]); % KMD ratios  correspond to CH2,H2,COO,CH2O,C2H4O,O,H2O,NH3,HCl series
    Filter_KMD_Series_Threshold = 3; % Minimum number of points on a KMD series to be used for refinement
    Filter_KMD_Unlimited_Iterations = true; % Enable or disable for unlimited iterations (true or false)   

    % Stage 5: Composition filter
    Filter_Composition = true; % Enable or disable filter (true or false)

    % Stage 6: Error filter
    Filter_Error = true; % Enable or disable filter (true or false)
    end
end