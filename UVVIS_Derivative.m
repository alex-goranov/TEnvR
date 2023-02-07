function varargout = UVVIS_Derivative(filename_sample)
%% Description

% This code takes a UV-VIS spectrum (wavelength in column one, absorbance values in column two)
% and take its derivative. Higher order derivatives can be taken by applying the code multiple times. 
% The spectrum must be located in a .csv file (comma-separated values). 

% Examples:
% UVVIS_Derivative('Sample 3_Final.csv')
% UVVIS_Derivative('Sample 1_Final_deriv.csv')
% UVVIS_Derivative('Sample 1_Final_deriv_deriv.csv')

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

% Import
Data=dlmread(filename_sample); 
Wavelengths=Data(:,1);
Absorbances=Data(:,2);

%Take the derivative
dydx = gradient(Absorbances(:)) ./ gradient(Wavelengths(:));


if nargout == 0
    dlmwrite([filename_sample(1:end-4) '_deriv.csv'],[Wavelengths,dydx])
    disp(['Finished applying UVVIS_Derivative on ' char(filename_sample)])
elseif nargout == 1
    disp(['Finished applying UVVIS_Derivative on ' char(filename_sample)])
    varargout{1}=[Wavelengths,dydx];
end

end