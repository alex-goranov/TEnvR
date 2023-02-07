function EEM_WaterRamanAverage(filename_FirstRep,range)
%% Description

% This code takes five fluorescence spectra (wavelength in column one, intensity values in column two)
% and averages them into one spectrum that is exported in a .csv file (comma-separated values). 

% This code is primarily used for averaging water Raman spectra (Ex=350 nm)
% but can be successfully used for any two-dimensional fluorescence spectra.

% Example:
% EEM_WaterRamanAverage('WaterRaman_1PM_1.txt','A38..B208')

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

clearvars -except filename_FirstRep range 
close all

filename_Rep1=filename_FirstRep;
filename_Rep2=[filename_FirstRep(1:end-5) '2.txt'];
filename_Rep3=[filename_FirstRep(1:end-5) '3.txt'];
filename_Rep4=[filename_FirstRep(1:end-5) '4.txt'];
filename_Rep5=[filename_FirstRep(1:end-5) '5.txt'];

Data1=readmatrix(filename_Rep1,'Range',range);
Data2=readmatrix(filename_Rep2,'Range',range);
Data3=readmatrix(filename_Rep3,'Range',range);
Data4=readmatrix(filename_Rep4,'Range',range);
Data5=readmatrix(filename_Rep5,'Range',range);


%% Stage 1: Data testing

if size(Data1,1) ~= size(Data2,1)
    error('Error! 1st and 2nd replicate are with different dimensions!')
end    
if size(Data1,1) ~= size(Data3,1)
    error('Error! 1st and 3nd replicate are with different dimensions!')
end    
if size(Data1,1) ~= size(Data4,1)
    error('Error! 1st and 4nd replicate are with different dimensions!')
end    
if size(Data1,1) ~= size(Data5,1)
    error('Error! 1st and 5nd replicate are with different dimensions!')
end  
    
%% Stage 2: Data averaging

Wavelengths=Data1(:,1);

if size(Wavelengths,1) ~= size(unique(Wavelengths),1)
   Wavelengths=transpose(Wavelengths(1):((Wavelengths(end)-Wavelengths(1))/(size(Wavelengths,1)-1)):Wavelengths(end));
end

Intensities=[Data1(:,2),Data2(:,2),Data3(:,2),Data4(:,2),Data5(:,2)];
AverageInt=mean(Intensities,2);

FinalMatrix=[Wavelengths,AverageInt];
dlmwrite([filename_FirstRep(1:end-6) '_AVG.csv'],FinalMatrix);

disp(['Finished applying EEM_WaterRamanAverage on ' char(filename_FirstRep) ' to ' char(filename_Rep5) '.csv'])

disp(['New file ' [filename_FirstRep(1:end-6) '_AVG.csv'] ' created!'])
end