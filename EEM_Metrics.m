%% Description

% This code takes one or multiple processed EEM spectra and calculates a
% variety of commonly used metrics (total spectral intensity, fluorescence index, humification index, etc.)
% The EEM spectrum must be located in a .csv file (comma-separated values). 

% Examples:
% Run code as EEM_Metrics in Command Window. In the pop-up window, use:
% filename_ref_Final.csv.csv - for running individual samples
% ALL - for running a whole dataset

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

global function_RunOnAllSamples
filename=inputdlg('Enter filename or use ALL','Filename',1);
filename=string(filename);
warning('off','MATLAB:xlswrite:AddSheet')
if filename == 'ALL'
    function_RunOnAllSamples=true;
    Files=dir('*_Final.csv'); % Specify file extension
    Files=natsortfiles(Files);
    for z=1:size(Files,1)
       FileName=Files(z).name;
       [~,~,Matrix1,Matrix2,Matrix3]=EEM_MetricsCalculation(FileName);
       M(z,:)=[string(Matrix1),string(Matrix2),string(Matrix3)];
    end
    
    Labels1={Files.name};
    Labels1=transpose(Labels1);
    output=[Labels1,M];

    Labels2={' ' 'Before Interpolation' 'Before Interpolation' 'Before Interpolation' 'Before Interpolation' 'Before Interpolation' 'Before Interpolation' ...
    'After Interpolation' 'After Interpolation' 'After Interpolation' 'After Interpolation' 'After Interpolation' 'After Interpolation' ...
    'Index' 'Index' 'Index' 'Index' 'Index' 'Index' 'Index' 'Index' 'Index' 'Index' 'Index' 'Index' 'Index' ...
    'Peak Alpha' 'Peak Alpha' 'Peak Alpha' 'Peak Alpha' 'Peak Alpha' 'Peak Alpha' ...
    'Peak Beta' 'Peak Beta' 'Peak Beta' 'Peak Beta' 'Peak Beta' 'Peak Beta'...
    'Tyrosine-A, def1' 'Tyrosine-A, def1' 'Tyrosine-A, def1' 'Tyrosine-A, def1' 'Tyrosine-A, def1' 'Tyrosine-A, def1'...
    'Tyrosine-A, def2' 'Tyrosine-A, def2' 'Tyrosine-A, def2' 'Tyrosine-A, def2' 'Tyrosine-A, def2' 'Tyrosine-A, def2'...
    'Tyrosine-B, def1' 'Tyrosine-B, def1' 'Tyrosine-B, def1' 'Tyrosine-B, def1' 'Tyrosine-B, def1' 'Tyrosine-B, def1'...
    'Tyrosine-B, def2' 'Tyrosine-B, def2' 'Tyrosine-B, def2' 'Tyrosine-B, def2' 'Tyrosine-B, def2' 'Tyrosine-B, def2'...
    'Tryptophan-A, def1' 'Tryptophan-A, def1' 'Tryptophan-A, def1' 'Tryptophan-A, def1' 'Tryptophan-A, def1' 'Tryptophan-A, def1'...
    'Tryptophan-A, def2' 'Tryptophan-A, def2' 'Tryptophan-A, def2' 'Tryptophan-A, def2' 'Tryptophan-A, def2' 'Tryptophan-A, def2'...
    'Tryptophan-B, def1' 'Tryptophan-B, def1' 'Tryptophan-B, def1' 'Tryptophan-B, def1' 'Tryptophan-B, def1' 'Tryptophan-B, def1'...
    'Tryptophan-B, def2' 'Tryptophan-B, def2' 'Tryptophan-B, def2' 'Tryptophan-B, def2' 'Tryptophan-B, def2' 'Tryptophan-B, def2'...
    'Peak N, def1' 'Peak N, def1' 'Peak N, def1' 'Peak N, def1' 'Peak N, def1' 'Peak N, def1'...
    'Peak N, def2' 'Peak N, def2' 'Peak N, def2' 'Peak N, def2' 'Peak N, def2' 'Peak N, def2'...
    'Peak D, def1' 'Peak D, def1' 'Peak D, def1' 'Peak D, def1' 'Peak D, def1' 'Peak D, def1'...
    'Peak D, def2' 'Peak D, def2' 'Peak D, def2' 'Peak D, def2' 'Peak D, def2' 'Peak D, def2'...
    'Peak E, def1' 'Peak E, def1' 'Peak E, def1' 'Peak E, def1' 'Peak E, def1' 'Peak E, def1'...
    'Peak E, def2' 'Peak E, def2' 'Peak E, def2' 'Peak E, def2' 'Peak E, def2' 'Peak E, def2'...   
    'Peak P, def1' 'Peak P, def1' 'Peak P, def1' 'Peak P, def1' 'Peak P, def1' 'Peak P, def1'...
    'Peak P, def2' 'Peak P, def2' 'Peak P, def2' 'Peak P, def2' 'Peak P, def2' 'Peak P, def2'...
    'Peak H' 'Peak H' 'Peak H' 'Peak H' 'Peak H' 'Peak H'...
    'Peak H mod' 'Peak H mod' 'Peak H mod' 'Peak H mod' 'Peak H mod' 'Peak H mod'...
    'Peak M-a' 'Peak M-a' 'Peak M-a' 'Peak M-a' 'Peak M-a' 'Peak M-a'...
    'Peak M-a mod' 'Peak M-a mod' 'Peak M-a mod' 'Peak M-a mod' 'Peak M-a mod' 'Peak M-a mod'...
    'Peak M-b def1' 'Peak M-b def1' 'Peak M-b def1' 'Peak M-b def1' 'Peak M-b def1' 'Peak M-b def1'...
    'Peak M-b def2' 'Peak M-b def2' 'Peak M-b def2' 'Peak M-b def2' 'Peak M-b def2' 'Peak M-b def2'...
    'Peak C-a def1' 'Peak C-a def1' 'Peak C-a def1' 'Peak C-a def1' 'Peak C-a def1' 'Peak C-a def1'...
    'Peak C-a def2' 'Peak C-a def2' 'Peak C-a def2' 'Peak C-a def2' 'Peak C-a def2' 'Peak C-a def2'...
    'Peak C-b def1' 'Peak C-b def1' 'Peak C-b def1' 'Peak C-b def1' 'Peak C-b def1' 'Peak C-b def1'...
    'Peak C-b def2' 'Peak C-b def2' 'Peak C-b def2' 'Peak C-b def2' 'Peak C-b def2' 'Peak C-b def2'...
    'Peak C+a' 'Peak C+a' 'Peak C+a' 'Peak C+a' 'Peak C+a' 'Peak C+a'...
    'Peak C+a mod' 'Peak C+a mod' 'Peak C+a mod' 'Peak C+a mod' 'Peak C+a mod' 'Peak C+a mod'...
    'Peak C+b' 'Peak C+b' 'Peak C+b' 'Peak C+b' 'Peak C+b' 'Peak C+b'};

    Labels3= {' ' 'Total Int' 'Average Int' 'Maximum Int' 'Peak (Ex)' 'Peak (Em)' 'Int around Peak' ...
    'Total Int' 'Average Int' 'Maximum Int' 'Peak (Ex)' 'Peak (Em)' 'Int around Peak'...
    'HIX_SYN_fulvics1' 'HIX_SYN_fulvics2' 'HIX_SYN_DOM' 'HIX_EM_old' 'HIX_EM_new' 'BIX_old' ...
    'BIX_new' 'FI_old' 'FI_new' 'T/C version1' 'T/C version2' 'Proctor Index' 'Perrette Index'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'...
    'Maximum Int' 'Sum of Int' 'Int Avg' 'Int Avg/Avg' 'Peak (Ex)' 'Peak (Em)'};

    output=[Labels2;Labels3;output];
    xlswrite('EEM Metrics Report.xlsx',output,1)
    disp('EEM_Metrics - Complete!')
else
    function_RunOnAllSamples=false;
    [Data,Data_Interp,Matrix1,Matrix2,Matrix3]=EEM_MetricsCalculation(filename);
    filename=char(filename);
    filename_new=filename(1:end-4);
    xlswrite([filename_new '_Metrics.xlsx'],Data,1,'A1');
    xlswrite([filename_new '_Metrics.xlsx'],Data_Interp,'Interpolated Data','A1');
    
    Labels1 = {' ' 'Before Interpolation' 'After Interpolation'};
    Labels2 = {'Total Spectral Intensity';'Average Spectral Intensity';'Maximum Spectral Intensity';...
        'Most Intense Peak (Ex)'; 'Most Intense Peak (Em)'; 'Sum of Intensity of Around Most Intense Peak'};
    output1=[Labels1;Labels2,string(Matrix1)];
    
    Labels3 ={'Index';'HIX_SYN_a'; 'HIX_SYN_b'; 'HIX_SYN_c'; 'HIX_EM_a'; 'HIX_EM_b'; 'BIX_old'; 'BIX_new'; 'FI_old'; 'FI_new'; 'T/C version a'; 'T/C version b'; 'Proctor Index'; 'Perrette Index'};
    output2=['Value';string(Matrix2)];
    output2=[Labels3,output2];
    
    Labels4 = {' ' 'Maximum Intensity' 'Summated Intensity' 'Averaged Intensity' 'Int Avg/Avg' 'Peak Position of Maximum  (Ex)' 'Peak Position of Maximum  (Em)'};
    Labels5 = {'Peak Alpha';'Peak Beta';'Tyrosine-1, def1';'Tyrosine-1, def2';'Tyrosine-2, def1';'Tyrosine-2, def2';'Tryptophan-1, def1';...
    'Tryptophan-1, def2';'Tryptophan-2, def1';'Tryptophan-2, def2';'Peak N, def1';'Peak N, def2';'Peak D, def1';'Peak D, def2';'Peak E, def1';'Peak E, def2';...
    'Peak P, def1';'Peak P, def2';'Peak H, def1';'Peak H, def2'; 'Peak M-1, def1';'Peak M-1, def2';'Peak M-2, def1'; 'Peak M-2, def2';'Peak C-1, def1';'Peak C-1, def2';...
    'Peak C-2, def1';'Peak C-2, def2';'Peak C(+)-1, def1'; 'Peak C(+)-1, def2'; 'Peak C(+)-2'};
    output3 = [Labels4; Labels5, string(Matrix3)];

    xlswrite([filename_new '_Metrics.xlsx'],output1,'Metrics','A1');
    xlswrite([filename_new '_Metrics.xlsx'],output2,'Metrics','F1');
    xlswrite([filename_new '_Metrics.xlsx'],output3,'Metrics','I1');
end

%% Internal functions 

% Function for calculating metrics on one EEM spectrum
function [Data,Data_Interp,Matrix1,Matrix2,Matrix3]=EEM_MetricsCalculation(filename)
global function_RunOnAllSamples

%% Stage 1: Import and organize data

Data=csvread(filename); % Import data
Ex=Data(1,:); Em=Data(:,1); % Extract excitation and emission ranges
Ex=Ex(2:end); Em=Em(2:end); % Extract excitation and emission ranges
Int=Data; Int(1,:)=[]; Int(:,1)=[]; % Extract intensity data 
Int=fillmissing(Int,'constant',0); % Any NaN values are turned into 0

% Make sure the dimensions are consistent
Size_Em=size(Em);
Size_Ex=size(Ex);
Size_Interp=size(Int);

if Size_Em(1)~= Size_Interp(1)
    error('Error! Check the dimensions of the emission of sample 1')
end
if Size_Ex(2)~= Size_Interp(2)
    error('Error! Check the dimensions of the excitation of sample 1')
end

%% Stage 2: Interpolate using interp2(x,y,v,xq,yq)

% Determine if interpolation is necessary

if Ex(2)-Ex(1) ~=1 | Em(2)-Em(1) ~=1
    % Creates X and Y matrices of the old EX and EM wavelengths
    [Interp_X,Interp_Y] = meshgrid(Ex,Em);

    % V is the old Int values
    Interp_V=Int;

    % Determine new Ex and EM dimensions
    Ex_Interp=round(Ex(1,1)):1:round(Ex(1,end));
    Em_Interp=round(Em(1,1)):1:round(Em(end,1));
    Em_Interp=transpose(Em_Interp);

    % Creates new X and Y matrices of the old EX and EM wavelengths
    [Interp_Xq,Interp_Yq] = meshgrid(Ex_Interp,Em_Interp);

    % Interpoloate
    Int_Interp=interp2(Interp_X,Interp_Y,Interp_V,Interp_Xq,Interp_Yq);
    Int_Interp=fillmissing(Int_Interp,'constant',0); %Any NaN values are turned into 0

    % Assemble new dataset
    Data_Interp=[0,Ex_Interp;Em_Interp,Int_Interp];
else
    Ex_Interp=Ex;
    Em_Interp=Em;
    Int_Interp=Int;
    Data_Interp=Data;    
end

% Check the data to make sure there are no negative intensity values
Variable_Negative_Elements_Interp=Int<0;
if sum(sum(Variable_Negative_Elements_Interp)) ~= 0
    error('Error! You have negative intensity values!')
end

Variable_Negative_Elements_Interp_Interp=Int_Interp<0;
if sum(sum(Variable_Negative_Elements_Interp_Interp)) ~= 0
    error('Error! You have negative intensity values!')
end

%% Stage 3: Calculate Metrics

Variable_Positive_Elements_Interp=Int>0; % Any NaN values are turned into 0
Variable_Positive_Elements_Interp_Interp=Int_Interp>0; % Any NaN values are turned into 0

% METRIC_1 = Total Spectral Intensity before data interpolation 
METRIC1_TotalIntensity=sum(sum(Int));

% METRIC_2 = Total Spectral Intensity after data interpolation 
METRIC2_TotalIntensity_Interp=sum(sum(Int_Interp));

% METRIC_3 = Average Spectral Intensity before data interpolation 
METRIC3_AverageIntensity=(sum(sum(Int)))/(sum(sum(Variable_Positive_Elements_Interp)));

% METRIC_4 = Average Spectral Intensity after data interpolation 
METRIC4_AverageIntensity_Interp=(sum(sum(Int_Interp)))/(sum(sum(Variable_Positive_Elements_Interp_Interp)));

% METRIC_5 = Maximum Spectral Intensity point before data interpolation
METRIC5_SpectralMaximumIntensity=max(max(Int));

% METRIC_6 = Maximum Spectral Intensity point after data interpolation 
METRIC6_SpectralMaximumIntensity_Interp=max(max(Int_Interp));

% METRIC 7a,b = Peak Position(Ex & Em) of Maximum Spectral Intensity before data interpolation 
[~,~,~,METRIC7a_SpectralMaximumPositionEX,METRIC7b_SpectralMaximumPositionEM] = IntensityRangeFunction(min(Ex),max(Ex),min(Em),max(Em),Data,Ex,Em);

% METRIC 8 = Sum of Intensity around (+-5) maximum peak position before data interpolation 
[~,METRIC8_SpectralMaximumIntensityRegion,~] = IntensityRangeFunction(METRIC7a_SpectralMaximumPositionEX-5,...
    METRIC7a_SpectralMaximumPositionEX+5,METRIC7b_SpectralMaximumPositionEM-5,METRIC7b_SpectralMaximumPositionEM+5,Data,Ex,Em);

% METRIC 9a,b = Peak Position(Ex & Em) of Maximum Spectral Intensity after data interpolation
[~,~,~,METRIC9a_SpectralMaximumPositionEX_Interp,METRIC9b_SpectralMaximumPositionEM_Interp] = IntensityRangeFunction(min(Ex_Interp),max(Ex_Interp),min(Em_Interp),max(Em_Interp),Data_Interp,Ex_Interp,Em_Interp);

% METRIC 10 = Sum of Intensity around (+-5) maximum peak position after data interpolation
[~,METRIC10_SpectralMaximumIntensityRegion_Interp,~] = IntensityRangeFunction(METRIC9a_SpectralMaximumPositionEX_Interp-5,...
    METRIC9a_SpectralMaximumPositionEX_Interp+5,METRIC9b_SpectralMaximumPositionEM_Interp-5,METRIC9b_SpectralMaximumPositionEM_Interp+5,Data_Interp,Ex_Interp,Em_Interp);

% METRIC 11a,b,c = Humification Index based on synchronous scans.
% a and b are for fulvic acids, c is for whole water samples. 

METRIC11a_HIXsyn=IntensitySearch_OnePoint(Data_Interp,470,488)/IntensitySearch_OnePoint(Data_Interp,360,378);
METRIC11b_HIXsyn=IntensitySearch_OnePoint(Data_Interp,400,418)/IntensitySearch_OnePoint(Data_Interp,360,378);
METRIC11c_HIXsyn=IntensitySearch_OnePoint(Data_Interp,390,408)/IntensitySearch_OnePoint(Data_Interp,355,373);

% METRIC 12a,b = Humification Index based on regular fluorescence
% aquisition. a is original, b is modified by Ohno. 
METRIC12a_HIXem=AreaUnderCurve(254,435,480,Data_Interp)/AreaUnderCurve(254,300,345,Data_Interp);
METRIC12b_HIXem=AreaUnderCurve(254,435,480,Data_Interp)/(AreaUnderCurve(254,300,345,Data_Interp)+AreaUnderCurve(254,435,480,Data_Interp));

% METRIC 13a,b,c,d,e = Peak Alpha. a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC13a_PeakAlpha_Max,METRIC13b_PeakAlpha_Sum,METRIC13c_PeakAlpha_Avg,METRIC13d_PeakAlpha_MaxEx,...
    METRIC13e_PeakAlpha_MaxEm]=IntensityRangeFunction(330,350,420,480,Data_Interp,Ex_Interp,Em_Interp);

M13=[METRIC13a_PeakAlpha_Max,METRIC13b_PeakAlpha_Sum,METRIC13c_PeakAlpha_Avg,METRIC13c_PeakAlpha_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC13d_PeakAlpha_MaxEx,METRIC13e_PeakAlpha_MaxEm];

% METRIC 14a,b,c,d,e = Peak Beta. a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC14a_PeakBeta_Max,METRIC14b_PeakBeta_Sum,METRIC14c_PeakBeta_Avg,METRIC14d_PeakBeta_MaxEx,METRIC14e_PeakBeta_MaxEm,...
    ]=IntensityRangeFunction(310,320,380,420,Data_Interp,Ex_Interp,Em_Interp);

M14=[METRIC14a_PeakBeta_Max,METRIC14b_PeakBeta_Sum,METRIC14c_PeakBeta_Avg,METRIC14c_PeakBeta_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC14d_PeakBeta_MaxEx,METRIC14e_PeakBeta_MaxEm];

% METRIC 15 = Freshness Index (BIX), original (peak beta/peak alpha)
METRIC15_BIXold=METRIC14a_PeakBeta_Max/METRIC13a_PeakAlpha_Max;

% METRIC 16 = Freshness Index (BIX), modified
[Variable_METRIC16Denominator,~,~,~,~] = IntensityRangeFunction(310,310,420,435,Data_Interp,Ex_Interp,Em_Interp);
METRIC16_BIXnew=IntensitySearch_OnePoint(Data_Interp,310,380)/Variable_METRIC16Denominator;

% METRIC 17 = Fluorescence Index, original definition as ratio of emission at
% 450 nm and 500 nm at EX=370 nm. For non-insturument-corrected spectra.
METRIC17_FIold=IntensitySearch_OnePoint(Data_Interp,370,450)/IntensitySearch_OnePoint(Data_Interp,370,500);

% METRIC 18 = Fluorescence Index, modified definition as ratio of emission at
% 470 nm and 520 nm at EX=370 nm. For insturument-corrected spectra.
METRIC18_FInew=IntensitySearch_OnePoint(Data_Interp,370,470)/IntensitySearch_OnePoint(Data_Interp,370,520);

% METRIC 19a,b = PeakT/PeakC ratio
% a is based on the original definition, where Peak T is emission at ex275/em350
% b is based on the definition of Peak T and Peak Cfrom Table 9.2
[Variable_METRIC19aDenominator,~,~,~,~] = IntensityRangeFunction(320,340,410,430,Data_Interp,Ex_Interp,Em_Interp);

METRIC19a_TtoC=IntensitySearch_OnePoint(Data_Interp,275,350)/Variable_METRIC19aDenominator;

[Variable_METRIC19bNominator,~,~,~,~] = IntensityRangeFunction(270,280,320,350,Data_Interp,Ex_Interp,Em_Interp);
[Variable_METRIC19bDenominator,~,~,~,~] = IntensityRangeFunction(330,350,420,480,Data_Interp,Ex_Interp,Em_Interp);

METRIC19b_TtoC=Variable_METRIC19bNominator/Variable_METRIC19bDenominator;

% METRIC 20 = Proctor Index, developed by Proctor et al. 2000

METRIC20_ProctorIndex=IntensitySearch_OnePoint(Data_Interp,350,420)/IntensitySearch_OnePoint(Data_Interp,390,470);

% METRIC 21 = Perrette Index, developed by Perrette et al. 2005
METRIC21_PerretteIndex=IntensitySearch_OnePoint(Data_Interp,364,514)/IntensitySearch_OnePoint(Data_Interp,364,457);

% METRIC 22 = Peak Tyrosine-1, original definition from Table 3.1
METRIC22_PeakTyrosine1=IntensitySearch_OnePoint(Data_Interp,230,305);

M22=[METRIC22_PeakTyrosine1,NaN,NaN,NaN,NaN,NaN];

% METRIC 23a,b,c,d,e = Peak Tyrosine-1, modified definition by Alex Goranov to include a range(+-5 nm)
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC23a_PeakTyrosine1_Max,METRIC23b_PeakTyrosine1_Sum,METRIC23c_PeakTyrosine1_Avg,METRIC23d_PeakTyrosine1_MaxEx,...
    METRIC23e_PeakTyrosine1_MaxEm] = IntensityRangeFunction(230-5,230+5,305-5,305+5,Data_Interp,Ex_Interp,Em_Interp);

M23=[METRIC23a_PeakTyrosine1_Max,METRIC23b_PeakTyrosine1_Sum,METRIC23c_PeakTyrosine1_Avg,...
    METRIC23c_PeakTyrosine1_Avg/METRIC4_AverageIntensity_Interp,METRIC23d_PeakTyrosine1_MaxEx,METRIC23e_PeakTyrosine1_MaxEm];

% METRIC 24 = Peak Tyrosine-2 (Peak B), original definition from Table 3.1
METRIC24_Tyrosine2=IntensitySearch_OnePoint(Data_Interp,275,305);

M24=[METRIC24_Tyrosine2,NaN,NaN,NaN,NaN,NaN];

% METRIC 25a,b,c,d,e = Peak Tyrosine-2 (Peak B), second original definition from Table 9.2
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)

[METRIC25a_PeakTyrosine2_Max,METRIC25b_PeakTyrosine2_Sum,METRIC25c_PeakTyrosine2_Avg,METRIC25d_PeakTyrosine2_MaxEx,...
    METRIC25e_PeakTyrosine2_MaxEm] = IntensityRangeFunction(270,280,300,320,Data_Interp,Ex_Interp,Em_Interp);

M25=[METRIC25a_PeakTyrosine2_Max,METRIC25b_PeakTyrosine2_Sum,METRIC25c_PeakTyrosine2_Avg,...
    METRIC25c_PeakTyrosine2_Avg/METRIC4_AverageIntensity_Interp,METRIC25d_PeakTyrosine2_MaxEx,METRIC25e_PeakTyrosine2_MaxEm];

% METRIC 26 = Peak Tryptophan-1, original definition from Table 3.1
METRIC26_PeakTryptophan1=IntensitySearch_OnePoint(Data_Interp,230,340);

M26=[METRIC26_PeakTryptophan1,NaN,NaN,NaN,NaN,NaN];

% METRIC 27a,b,c,d,e = Peak Tryptophan-1, modified definition by Alex Goranov to include a range(+-5 nm)
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC27a_PeakTryptophan1_Max,METRIC27b_PeakTryptophan1_Sum,METRIC27c_PeakTryptophan1_Avg,METRIC27d_PeakTryptophan1_MaxEx,...
    METRIC27e_PeakTryptophan1_MaxEm] = IntensityRangeFunction(230-5,230+5,340-5,340+5,Data_Interp,Ex_Interp,Em_Interp);

M27=[METRIC27a_PeakTryptophan1_Max,METRIC27b_PeakTryptophan1_Sum,METRIC27c_PeakTryptophan1_Avg,...
    METRIC27c_PeakTryptophan1_Avg/METRIC4_AverageIntensity_Interp,METRIC27d_PeakTryptophan1_MaxEx,METRIC27e_PeakTryptophan1_MaxEm];

% METRIC 28 = Peak Tryptophan-2, original definition from Table 3.1
METRIC28_PeakTryptophan2=IntensitySearch_OnePoint(Data_Interp,275,340);

M28=[METRIC28_PeakTryptophan2,NaN,NaN,NaN,NaN,NaN];

% METRIC 29a,b,c,d,e = Peak Tryptophan-2 (Peak T), second original definition from Table 9.2
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC29a_PeakTryptophan2_Max,METRIC29b_PeakTryptophan2_Sum,METRIC29c_PeakTryptophan2_Avg,METRIC29d_PeakTryptophan2_MaxEx,...
    METRIC29e_PeakTryptophan2_MaxEm] = IntensityRangeFunction(270,280,320,350,Data_Interp,Ex_Interp,Em_Interp);

M29=[METRIC29a_PeakTryptophan2_Max,METRIC29b_PeakTryptophan2_Sum,METRIC29c_PeakTryptophan2_Avg,...
    METRIC29c_PeakTryptophan2_Avg/METRIC4_AverageIntensity_Interp,METRIC29d_PeakTryptophan2_MaxEx,METRIC29e_PeakTryptophan2_MaxEm];

% METRIC 30 = Peak N, original definition from Table 3.1
METRIC30_PeakN=IntensitySearch_OnePoint(Data_Interp,280,370);

M30=[METRIC30_PeakN,NaN,NaN,NaN,NaN,NaN];

% METRIC 31a,b,c,d,e = Peak N, modified definition by Alex Goranov to include a range(+-5 nm)
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC31a_PeakN_Max,METRIC31b_PeakN_Sum,METRIC31c_PeakN_Avg,METRIC31d_PeakN_MaxEx,...
    METRIC31e_PeakN_MaxEm] = IntensityRangeFunction(280-5,280+5,370-5,370+5,Data_Interp,Ex_Interp,Em_Interp);

M31=[METRIC31a_PeakN_Max,METRIC31b_PeakN_Sum,METRIC31c_PeakN_Avg,METRIC31c_PeakN_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC31d_PeakN_MaxEx,METRIC31e_PeakN_MaxEm];

% METRIC 32 = Peak D, original definition from Table 2.1
METRIC32_PeakD=IntensitySearch_OnePoint(Data_Interp,390,509);

M32=[METRIC32_PeakD,NaN,NaN,NaN,NaN,NaN];

% METRIC 33a,b,c,d,e = Peak D, modified definition by Alex Goranov to include a range(+-5 nm)
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC33a_PeakD_Max,METRIC33b_PeakD_Sum,METRIC33c_PeakD_Avg,METRIC33d_PeakD_MaxEx,...
    METRIC33e_PeakD_MaxEm] = IntensityRangeFunction(390-5,390+5,509-5,509+5,Data_Interp,Ex_Interp,Em_Interp);

M33=[METRIC33a_PeakD_Max,METRIC33b_PeakD_Sum,METRIC33c_PeakD_Avg,METRIC33c_PeakD_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC33d_PeakD_MaxEx,METRIC33e_PeakD_MaxEm];

% METRIC 34 = Peak E, original definition from Table 2.1
METRIC34_PeakE=IntensitySearch_OnePoint(Data_Interp,455,521);

M34=[METRIC34_PeakE,NaN,NaN,NaN,NaN,NaN];

% METRIC 35a,b,c,d,e = Peak E, modified definition by Alex Goranov to include a range(+-5 nm)
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC35a_PeakE_Max,METRIC35b_PeakE_Sum,METRIC35c_PeakE_Avg,METRIC35d_PeakE_MaxEx,...
    METRIC35e_PeakE_MaxEm] = IntensityRangeFunction(455-5,455+5,521-5,521+5,Data_Interp,Ex_Interp,Em_Interp);

M35=[METRIC35a_PeakE_Max,METRIC35b_PeakE_Sum,METRIC35c_PeakE_Avg,METRIC35c_PeakE_Avg/METRIC4_AverageIntensity_Interp...
    METRIC35d_PeakE_MaxEx,METRIC35e_PeakE_MaxEm];

% METRIC 36 = Peak P, Pigment-like, original definition from Table 3.1
METRIC36_PeakP=IntensitySearch_OnePoint(Data_Interp,398,660);

M36=[METRIC36_PeakP,NaN,NaN,NaN,NaN,NaN];

% METRIC 37a,b,c,d,e = Peak P, modified definition by Alex Goranov to include a range(+-5 nm)
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC37a_PeakP_Max,METRIC37b_PeakP_Sum,METRIC37c_PeakP_Avg,METRIC37d_PeakP_MaxEx,...
    METRIC37e_PeakP_MaxEm] = IntensityRangeFunction(398-5,398+5,660-5,660+5,Data_Interp,Ex_Interp,Em_Interp);

M37=[METRIC37a_PeakP_Max,METRIC37b_PeakP_Sum,METRIC37c_PeakP_Avg,METRIC37c_PeakP_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC37d_PeakP_MaxEx,METRIC37e_PeakP_MaxEm];

% METRIC 38a,b,c,d,e = Peak H, Photobleached, original definition from Table 3.1
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC38a_PeakHa_Max,METRIC38b_PeakHa_Sum,METRIC38c_PeakHa_Avg,METRIC38d_PeakHa_MaxEx,...
    METRIC38e_PeakHa_MaxEm] = IntensityRangeFunction(230,230,275,350,Data_Interp,Ex_Interp,Em_Interp);

M38=[METRIC38a_PeakHa_Max,METRIC38b_PeakHa_Sum,METRIC38c_PeakHa_Avg,METRIC38c_PeakHa_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC38d_PeakHa_MaxEx,METRIC38e_PeakHa_MaxEm];

% METRIC 39a,b,c,d,e = Peak H, Photobleached, modified definition by Alex Goranov to include a range(+-5 nm)
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC39a_PeakHb_Max,METRIC39b_PeakHb_Sum,METRIC39c_PeakHb_Avg,METRIC39d_PeakHb_MaxEx,...
    METRIC39e_PeakHb_MaxEm] = IntensityRangeFunction(230-5,230+5,275,350,Data_Interp,Ex_Interp,Em_Interp);

M39=[METRIC39a_PeakHb_Max,METRIC39b_PeakHb_Sum,METRIC39c_PeakHb_Avg,METRIC39c_PeakHb_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC39d_PeakHb_MaxEx,METRIC39e_PeakHb_MaxEm];

% METRIC 40a,b,c,d,e = Peak M-1, original definition from Table 3.1
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC40a_PeakM1a_Max,METRIC40b_PeakM1a_Sum,METRIC40c_PeakM1a_Avg,METRIC40d_PeakM1a_MaxEx,...
    METRIC40e_PeakM1a_MaxEm] = IntensityRangeFunction(240,240,350,400,Data_Interp,Ex_Interp,Em_Interp);

M40=[METRIC40a_PeakM1a_Max,METRIC40b_PeakM1a_Sum,METRIC40c_PeakM1a_Avg,METRIC40c_PeakM1a_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC40d_PeakM1a_MaxEx,METRIC40e_PeakM1a_MaxEm];

% METRIC 41a,b,c,d,e = Peak M-1, modified definition by Alex Goranov to include a range(+-5 nm)
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC41a_PeakM1b_Max,METRIC41b_PeakM1b_Sum,METRIC41c_PeakM1b_Avg,METRIC41d_PeakM1b_MaxEx,...
    METRIC41e_PeakM1b_MaxEm] = IntensityRangeFunction(240-5,240+5,350,400,Data_Interp,Ex_Interp,Em_Interp);

M41=[METRIC41a_PeakM1b_Max,METRIC41b_PeakM1b_Sum,METRIC41c_PeakM1b_Avg,METRIC41c_PeakM1b_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC41d_PeakM1b_MaxEx,METRIC41e_PeakM1b_MaxEm];

% METRIC 42a,b,c,d,e = Peak M-2 (Peak M), original definition from Table 3.1
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC42a_PeakM2a_Max,METRIC42b_PeakM2a_Sum,METRIC42c_PeakM2a_Avg,METRIC42d_PeakM2a_MaxEx,...
    METRIC42e_PeakM2a_MaxEm] = IntensityRangeFunction(290,310,370,420,Data_Interp,Ex_Interp,Em_Interp);

M42=[METRIC42a_PeakM2a_Max,METRIC42b_PeakM2a_Sum,METRIC42c_PeakM2a_Avg,METRIC42c_PeakM2a_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC42d_PeakM2a_MaxEx,METRIC42e_PeakM2a_MaxEm];

% METRIC 43a,b,c,d,e = Peak M-2 (Peak M), second original definition from Table 9.2
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC43a_PeakM2b_Max,METRIC43b_PeakM2b_Sum,METRIC43c_PeakM2b_Avg,METRIC43d_PeakM2b_MaxEx,...
    METRIC43e_PeakM2b_MaxEm] = IntensityRangeFunction(310,320,380,420,Data_Interp,Ex_Interp,Em_Interp);

M43=[METRIC43a_PeakM2b_Max,METRIC43b_PeakM2b_Sum,METRIC43c_PeakM2b_Avg,METRIC43c_PeakM2b_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC43d_PeakM2b_MaxEx,METRIC43e_PeakM2b_MaxEm];

% METRIC 44a,b,c,d,e = Peak C-1 (Peak A), original definition from Table 3.1
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC44a_PeakC1a_Max,METRIC44b_PeakC1a_Sum,METRIC44c_PeakC1a_Avg,METRIC44d_PeakC1a_MaxEx,...
    METRIC44e_PeakC1a_MaxEm] = IntensityRangeFunction(260,260,400,460,Data_Interp,Ex_Interp,Em_Interp);

M44=[METRIC44a_PeakC1a_Max,METRIC44b_PeakC1a_Sum,METRIC44c_PeakC1a_Avg,METRIC44c_PeakC1a_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC44d_PeakC1a_MaxEx,METRIC44e_PeakC1a_MaxEm];

% METRIC 45a,b,c,d,e = Peak C-1 (Peak A), second original definition from Table 9.2
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC45a_PeakC1b_Max,METRIC45b_PeakC1b_Sum,METRIC45c_PeakC1b_Avg,METRIC45d_PeakC1b_MaxEx,...
    METRIC45e_PeakC1b_MaxEm] = IntensityRangeFunction(250,260,380,480,Data_Interp,Ex_Interp,Em_Interp);

M45=[METRIC45a_PeakC1b_Max,METRIC45b_PeakC1b_Sum,METRIC45c_PeakC1b_Avg,METRIC45c_PeakC1b_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC45d_PeakC1b_MaxEx,METRIC45e_PeakC1b_MaxEm];

% METRIC 46a,b,c,d,e = Peak C-2 (Peak C), original definition from Table 3.1
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC46a_PeakC2a_Max,METRIC46b_PeakC2a_Sum,METRIC46c_PeakC2a_Avg,METRIC46d_PeakC2a_MaxEx,...
    METRIC46e_PeakC2a_MaxEm] = IntensityRangeFunction(320,365,420,470,Data_Interp,Ex_Interp,Em_Interp);

M46=[METRIC46a_PeakC2a_Max,METRIC46b_PeakC2a_Sum,METRIC46c_PeakC2a_Avg,METRIC46c_PeakC2a_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC46d_PeakC2a_MaxEx,METRIC46e_PeakC2a_MaxEm];

% METRIC 47a,b,c,d,e = Peak C-2 (Peak C), second original definition from Table 9.2
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC47a_PeakC2b_Max,METRIC47b_PeakC2b_Sum,METRIC47c_PeakC2b_Avg,METRIC47d_PeakC2b_MaxEx,...
    METRIC47e_PeakC2b_MaxEm] = IntensityRangeFunction(330,350,420,480,Data_Interp,Ex_Interp,Em_Interp);

M47=[METRIC47a_PeakC2b_Max,METRIC47b_PeakC2b_Sum,METRIC47c_PeakC2b_Avg,METRIC47c_PeakC2b_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC47d_PeakC2b_MaxEx,METRIC47e_PeakC2b_MaxEm];

% METRIC 48a,b,c,d,e = Peak C(+)-1, original definition from Table 3.1
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC48a_PeakCplus1a_Max,METRIC48b_PeakCplus1a_Sum,METRIC48c_PeakCplus1a_Avg,METRIC48d_PeakCplus1a_MaxEx,...
    METRIC48e_PeakCplus1a_MaxEm] = IntensityRangeFunction(250,250,470,504,Data_Interp,Ex_Interp,Em_Interp);

M48=[METRIC48a_PeakCplus1a_Max,METRIC48b_PeakCplus1a_Sum,METRIC48c_PeakCplus1a_Avg,METRIC48c_PeakCplus1a_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC48d_PeakCplus1a_MaxEx,METRIC48e_PeakCplus1a_MaxEm];

% METRIC 49a,b,c,d,e = Peak C(+)-1, modified definition by Alex Goranov to include a range(+-5 nm)
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC49a_PeakCplus1b_Max,METRIC49b_PeakCplus1b_Sum,METRIC49c_PeakCplus1b_Avg,METRIC49d_PeakCplus1b_MaxEx,...
    METRIC49e_PeakCplus1b_MaxEm] = IntensityRangeFunction(250-5,250+5,470,504,Data_Interp,Ex_Interp,Em_Interp);

M49=[METRIC49a_PeakCplus1b_Max,METRIC49b_PeakCplus1b_Sum,METRIC49c_PeakCplus1b_Avg,METRIC49c_PeakCplus1b_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC49d_PeakCplus1b_MaxEx,METRIC49e_PeakCplus1b_MaxEm];

% METRIC 50a,b,c,d,e = Peak C(+)-2, original definition from Table 3.1
% a=maximum intensity, b=sum of intensities,c=average intensity
% d and e = peak position of the maximum intensity (Ex and Em, respectively)
[METRIC50a_PeakCplus2_Max,METRIC50b_PeakCplus2_Sum,METRIC50c_PeakCplus2_Avg,METRIC50d_PeakCplus2_MaxEx,...
    METRIC50e_PeakCplus2_MaxEm] = IntensityRangeFunction(385,420,470,504,Data_Interp,Ex_Interp,Em_Interp);

M50=[METRIC50a_PeakCplus2_Max,METRIC50b_PeakCplus2_Sum,METRIC50c_PeakCplus2_Avg,METRIC50c_PeakCplus2_Avg/METRIC4_AverageIntensity_Interp,...
    METRIC50d_PeakCplus2_MaxEx,METRIC50e_PeakCplus2_MaxEm];

%% Stage 4: Export Data

if function_RunOnAllSamples 
    Matrix1 = [METRIC1_TotalIntensity,METRIC3_AverageIntensity,METRIC5_SpectralMaximumIntensity,...
        METRIC7a_SpectralMaximumPositionEX,METRIC7b_SpectralMaximumPositionEM,METRIC8_SpectralMaximumIntensityRegion,...
        METRIC2_TotalIntensity_Interp,METRIC4_AverageIntensity_Interp,METRIC6_SpectralMaximumIntensity_Interp,...
        METRIC9a_SpectralMaximumPositionEX_Interp,METRIC9b_SpectralMaximumPositionEM_Interp,METRIC10_SpectralMaximumIntensityRegion_Interp];
    
    Matrix2 = [METRIC11a_HIXsyn;METRIC11b_HIXsyn;METRIC11c_HIXsyn;METRIC12a_HIXem;METRIC12b_HIXem;METRIC15_BIXold;METRIC16_BIXnew;METRIC17_FIold;METRIC18_FInew;METRIC19a_TtoC;METRIC19b_TtoC;METRIC20_ProctorIndex;METRIC21_PerretteIndex];
    Matrix2 = transpose(Matrix2);
    Matrix3 = [M13,M14,M22,M23,M24,M25,M26,M27,M28,M29,M30,M31,M32,M33,M34,M35,M36,M37,M38,M39,M40,...
        M41,M42,M43,M44,M45,M46,M47,M48,M49,M50];
else
    Matrix1 = [METRIC1_TotalIntensity,METRIC2_TotalIntensity_Interp;METRIC3_AverageIntensity,METRIC4_AverageIntensity_Interp;...
        METRIC5_SpectralMaximumIntensity,METRIC6_SpectralMaximumIntensity_Interp;METRIC7a_SpectralMaximumPositionEX,METRIC9a_SpectralMaximumPositionEX_Interp;...
        METRIC7b_SpectralMaximumPositionEM,METRIC9b_SpectralMaximumPositionEM_Interp;METRIC8_SpectralMaximumIntensityRegion,METRIC10_SpectralMaximumIntensityRegion_Interp];
    Matrix2 = [METRIC11a_HIXsyn;METRIC11b_HIXsyn;METRIC11c_HIXsyn;METRIC12a_HIXem;METRIC12b_HIXem;METRIC15_BIXold;METRIC16_BIXnew;METRIC17_FIold;METRIC18_FInew;METRIC19a_TtoC;METRIC19b_TtoC;METRIC20_ProctorIndex;METRIC21_PerretteIndex];
    Matrix3 = [M13;M14;M22;M23;M24;M25;M26;M27;M28;M29;M30;M31;M32;M33;M34;M35;M36;M37;M38;M39;M40;...
        M41;M42;M43;M44;M45;M46;M47;M48;M49;M50];
end

disp(['Finished applying EEM_Metrics on ' char(filename)])
end

% This function obtains the sum/averages/max of of intensities in a given range.
% Function also produces the position (Ex and EM of the highest intensity point in the given range)
function [Intensity_Max,Intensity_Sum,Intensity_Avg,METRIC_PositionEX,METRIC_PositionEM] = IntensityRangeFunction(EX_min,EX_max,EM_min,EM_max,Data,Ex,Em)

EX_range=Ex(Ex >= EX_min & Ex <= EX_max);
EM_range=Em(Em >= EM_min & Em <= EM_max);

if size(EM_range,1)~=0 && size(EX_range,2)~=0

for i=1:size(EX_range,2)
    for j=1:size(EM_range,1)
        A(j,i) = IntensitySearch_OnePoint(Data,EX_range(i),EM_range(j));
    end
end

Positive_Elements=A>0;
Intensity_Sum = sum(sum(A));
Intensity_Max = max(max(A));
Intensity_Avg = Intensity_Sum/sum(sum(Positive_Elements));

if Intensity_Max >0
    Variable_Row=find(any(Data==Intensity_Max,2));
    Variable_Column=find(any(Data==Intensity_Max,1));
    METRIC_PositionEX=Data(1,Variable_Column);
    METRIC_PositionEM=Data(Variable_Row,1);
    if length(METRIC_PositionEX) ~=1 % Make sure it is only 1 value in each direction
        METRIC_PositionEX=NaN;
        METRIC_PositionEM=NaN;
    end
    if length(METRIC_PositionEM) ~=1 % Make sure it is only 1 value in each direction
        METRIC_PositionEX=NaN;
        METRIC_PositionEM=NaN;
    end
else
    METRIC_PositionEX=NaN;
    METRIC_PositionEM=NaN;
end
else
    Intensity_Max=NaN;
    Intensity_Sum=NaN;
    Intensity_Avg=NaN;
    METRIC_PositionEX=NaN;
    METRIC_PositionEM=NaN;
end

end

% This function finds an intensity value in a "Data" matrix with given
% excitation and emission wavelengths
function Intensity_OnePoint = IntensitySearch_OnePoint(Data,Wavelength_Excitation,Wavelength_Emission)

[~,column]=ismember(Wavelength_Excitation,Data(1,:));
[~,row] = ismember(Wavelength_Emission,Data(:,1));

if column == 0
    Intensity_OnePoint = 0;
elseif row == 0
    Intensity_OnePoint = 0;
else
    Intensity_OnePoint = Data(row,column); 
end

end

% This function will find the area under the curve between EM_start and EM_end
% at a certain EX wavelength.
function Area = AreaUnderCurve(EX,EM_start,EM_end,Data)

EM_InterpegrationRange=EM_start:1:EM_end;

for j=1:size(EM_InterpegrationRange,2)
    B(j) = IntensitySearch_OnePoint(Data,EX,EM_InterpegrationRange(j));
end

Area=trapz(EM_InterpegrationRange,B);
end