function varargout = UVVIS_Metrics (filename)
%% Description

% This code takes a UV-VIS spectrum (wavelength in column one, absorbance values in column two)
% and calculates a variety of commonly used metrics (absorbance at 254, slope ratios, etc.)
% The spectrum must be located in a .csv file (comma-separated values). 

% Example:
% UVVIS_Metrics('Sample 1_Final.csv')

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

%% Data Import

close all
clearvars -except filename

Data_Decadic=dlmread(filename);
Data_Decadic=sortrows(Data_Decadic,1); % Sorts data short -> long wavelength 

%% Stage 1 = Calculate Metrics

Data_Napierian=[Data_Decadic(:,1),100*Data_Decadic(:,2)*(log(10))];
METRIC_Abs230Dec=Absorbance_OnePoint(Data_Decadic,230);
METRIC_Abs254Dec=Absorbance_OnePoint(Data_Decadic,254);
METRIC_Abs280Dec=Absorbance_OnePoint(Data_Decadic,280);
METRIC_alpha350=Absorbance_OnePoint(Data_Napierian,350);
METRIC_alpha325=Absorbance_OnePoint(Data_Napierian,325);
METRIC_alpha300=Absorbance_OnePoint(Data_Napierian,300);
METRIC_alpha250=Absorbance_OnePoint(Data_Napierian,250);
METRIC_alpha254=Absorbance_OnePoint(Data_Napierian,254);
METRIC_E4E6=Absorbance_OnePoint(Data_Decadic,465)/Absorbance_OnePoint(Data_Decadic,665);
METRIC_E2E3=Absorbance_OnePoint(Data_Decadic,250)/Absorbance_OnePoint(Data_Decadic,365);
METRIC_TotalCDOM=AreaUnderCurve(Data_Decadic,250,450);

[Abs_275295,Wav_275295]=Absorbance_Range(Data_Napierian,275,295);
[Abs_290320,Wav_290320]=Absorbance_Range(Data_Napierian,290,320);
[Abs_350400,Wav_350400]=Absorbance_Range(Data_Napierian,350,400);
[Abs_350450,Wav_350450]=Absorbance_Range(Data_Napierian,350,450);

[METRIC_Slope275295,~]=polyfit(Wav_275295,log(Abs_275295),1);
[METRIC_Slope290320,~]=polyfit(Wav_290320,log(Abs_290320),1);
[METRIC_Slope350400,~]=polyfit(Wav_350400,log(Abs_350400),1);
[METRIC_Slope350450,~]=polyfit(Wav_350450,log(Abs_350450),1);

METRIC_Slope275295=abs(METRIC_Slope275295(1));
METRIC_Slope290320=abs(METRIC_Slope290320(1));
METRIC_Slope350400=abs(METRIC_Slope350400(1));
METRIC_Slope350450=abs(METRIC_Slope350450(1));
METRIC_SlopeRatio=METRIC_Slope275295/METRIC_Slope350400;

%% Stage 2 = Metrics Export

matrix=[METRIC_Abs230Dec;METRIC_Abs254Dec; METRIC_Abs280Dec; METRIC_alpha350; METRIC_alpha325; METRIC_alpha300; METRIC_alpha250; METRIC_alpha254; METRIC_E4E6; METRIC_E2E3; METRIC_TotalCDOM; METRIC_Slope275295; METRIC_Slope290320;...
    METRIC_Slope350400;METRIC_Slope350450;METRIC_SlopeRatio];

if nargout == 0
    labels={'Abs-230 (dec)' 'Abs-254 (dec)' 'Abs-280 (dec)' 'Alpha350' 'Alpha325' 'Alpha300' 'Alpha250' 'Alpha254' 'E4toE6' 'E2toE3' 'Total CDOM' 'Slope(275-295)' 'Slope(290-320)' 'Slope(350-400)' 'Slope(350-450)' 'SR'};
    labels=transpose(labels);
    output=[labels,string(matrix)];
    warning('off','MATLAB:xlswrite:AddSheet');
    xlswrite([char(filename(1:end-4)) '_Metrics.xlsx'],{'Data Decadic'},1,'A1');
    xlswrite([char(filename(1:end-4)) '_Metrics.xlsx'],Data_Decadic,1,'A2');

    xlswrite([char(filename(1:end-4)) '_Metrics.xlsx'],{'Data Napierian'},1,'D1');    
    xlswrite([char(filename(1:end-4)) '_Metrics.xlsx'],Data_Napierian,1,'D2');

    xlswrite([char(filename(1:end-4)) '_Metrics.xlsx'],{char(filename(1:end-4))},1,'H1'); 
    xlswrite([char(filename(1:end-4)) '_Metrics.xlsx'],output,1,'G2'); 
elseif nargout == 3
    varargout{1}=matrix;
    varargout{2}=Data_Decadic;
    varargout{3}=Data_Napierian;
end
disp(['Finished applying UVVIS_Metrics on ' char(filename)])
end

%% Internal functions

% Extracts absorbance values over a range of wavelengths 
function [Abs_Range,Wav_Range] = Absorbance_Range(Data,Start,End)

Wav_Range=Start:1:End;

for j=1:size(Wav_Range,2)
    Abs_Range(j)=Absorbance_OnePoint(Data,Wav_Range(j));
end

end

% Extracts an absorbance value at a specified wavelength
function Abs_OnePoint = Absorbance_OnePoint(Data,Wavelength)

[~,row] = ismember(Wavelength,Data(:,1));
Abs_OnePoint = Data(row,2);
end

% Quantifies the area under a curve (integrates) data over a specified range
function Area = AreaUnderCurve(Data,Start,End)

IntegrationRange=Start:1:End;

for j=1:size(IntegrationRange,2)
    B(j) = Absorbance_OnePoint(Data,IntegrationRange(j));
end

Area=trapz(IntegrationRange,B);
end