%% Description

% This code takes multiple 1D NMR spectra (as exported by TopSpin)
% and processes them prior to PCA: 

    % 1) Reformats spectra into 2D arrays (chemical shift vs intensity)
    % 2) Denoises them (optional)
    % 3) Bins them (optional & tunable)
    % 4) Interpolates them (optional) 
    % 5) Aligns them into a matrix
    % 6) Normalizes them to total spectral intensity
    % 7) Exports everything in Excel. 

% Example: 
% Run the code by typing NMR_Automation in the command windwow.

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
    
%% Configuration
clear % Clears the workspace
close all % Closes any other MATLAB windows (e.g., Figure pop-ups)

% Define TopSpin version corresponding to file format and export command:
% topspin_new = false - convbin2ascii command is used, older TopSpin versions 
% topspin_new = true - totxt command is used, newer TopSpin versions 
topspin_new = true;

% Interpolate all data to make the same size
Interpolation=true;

% Spectral denoising - Identify signals that are 5*STD(noise) per Halouska et al., (2006) 
Denoise = false;

% Bin size = every binsize-number of points (e.g., every 200) will be averaged
Binning=true;
binsize=200;

%% Stage 1: Reformat files from .txt files from TopSpin into 2D arrays

Files=dir('*.txt'); % Specify file extension
for z=1:size(Files,1)
   FileName=Files(z).name;
   if topspin_new
       NMR_FormatNewTopspin(FileName);
   else
       NMR_FormatOldTopspin(FileName)
   end
end

%% Stage 2: Interpolation

Files=dir('*_ref.txt');
if Interpolation
    for z=1:size(Files,1)
       FileName=Files(z).name;
       Data=load(FileName);
       Data=sortrows(Data,1);    % Sort
       Data=Data(Data(:,1)>0,:); % Remove negative
       Endmember(z)=Data(end,1); % End point
       Sizes(z)=size(Data,1);
    end
    Max=max(Endmember);
    Size=mean(Sizes);
    
    for z=1:size(Files,1)
       FileName=Files(z).name;
       Data=load(FileName);
       
       Data=sortrows(Data,1);    % Sort
       Data=Data(Data(:,1)>0,:); % Remove negative

       [x,~,locMembers] = unique(Data(:,1),'rows','stable');
       v = splitapply(@(v) mean(v), Data(:,2), locMembers);
       
       xq=0:(Max/Size):Max;
       vq1 = interp1(x,v,xq,'linear');
       NewData=[xq',vq1'];
       
       dlmwrite([FileName(1:end-4) '_interp.txt'],NewData,'delimiter','\t','precision','%.3f');
       disp(['Finished interpolating ' char(FileName)])
    end
    Files=dir('*_interp.txt');
end

%% Stage 3: Denoising

if Denoise
    for z=1:size(Files,1)
       FileName=Files(z).name;
       NMR_Denoise(FileName);
    end
    Files=dir('*_denoised.txt');
end

%% Stage 4: Binning

if Binning
    for z=1:size(Files,1)
       FileName=Files(z).name;
       NMR_Bin(FileName,binsize);
    end
    Files=dir('*bin.txt');
end

%% Stage 5: Alignment & Normalization
    for z=1:size(Files,1)
       FileName=Files(z).name;
       Data=load(FileName);
       Data=sortrows(Data,1); %Sort
       Sizes(z)=size(Data,1); %number of points
    end

if std(Sizes) ~=0
    error('Error! The loaded spectra are of unequal size! Please enable the interpolation!')
end   

for z=1:size(Files,1)
    FileName=Files(z).name;
    Data=load(FileName);
    Data=sortrows(Data,1); %Sort
    Int=Data(:,2);
    M(z,:)=Int;
end
PPM=Data(:,1);
      
% replace negative and NaN values with zero
M(M<0)=0; 
M(isnan(M))=0;

% Renove PPM values with all zeros
indexes=find(sum(M ==0) == size(M,1));

M_original=M;
PPM_original=PPM;

M(:,indexes)=[];
PPM(indexes)=[];

% Normalize 
TotalSpectralInt=sum(M,2);
for i=1:size(M,1)
    RawInt=M(i,:);
    Norm(i,:)=RawInt./TotalSpectralInt(i);
end

     
%% Stage 6: Export

warning('off','MATLAB:xlswrite:AddSheet');

xlswrite('NMR Alignment Matrix.xlsx',{Files.name},1,'B1')
xlswrite('NMR Alignment Matrix.xlsx',PPM_original,1,'A2')
xlswrite('NMR Alignment Matrix.xlsx',M_original',1,'B2')

xlswrite('NMR Alignment Matrix.xlsx',{Files.name},'Normalized','B1')
xlswrite('NMR Alignment Matrix.xlsx',PPM,'Normalized','A2')
xlswrite('NMR Alignment Matrix.xlsx',Norm','Normalized','B2')


disp('NMR_Automation - Complete!')

%% Internal functions

% Function for converting new-version TopSpin exported files (using totxt) 
function NMR_FormatNewTopspin(filename)

C=importdata(filename);

% Export chemical shift axis parameters and figure out ppm dimensions
Size=C.textdata(6);
Size= str2double(regexp(Size,'\d+[\.]?\d*','match','once'));
B = strsplit(string(C.textdata(4)),'RIGHT');
Left=str2double(regexp(B(1),'[+-]?\d+[\.]?\d*','match','once'));
Right=str2double(regexp(B(2),'[+-]?\d+[\.]?\d*','match','once'));

% Recreate the chemical shift axis
PPM=[1:1:size(C.data,1)]';
PPM=PPM*((Left-Right)/Size);
PPM=PPM+Right;

% Consolidate intensity and axis arrays. Export.
Spectrum=[PPM,flip(C.data)];
dlmwrite([filename(1:end-4) '_ref.txt'],Spectrum);

disp(['Finished reformatting ' char(filename)])
end

% Function for converting old TopSpin exported files (using convbin2asc) 
function NMR_FormatOldTopspin(filename)
Data=load(filename);     

PPM=Data(:,4); % In old files, chemical shift data is in the 4th column        
Intensity=Data(:,2); % In old files, intensity data is in the 2nd column 

Output=[PPM, Intensity];
dlmwrite([filename(1:end-4) '_ref.txt'],Output,'\t');

disp(['Finished reformatting ' char(filename)])
end

% Function for denoising NMR spectra
function NMR_Denoise(filename)

% Identify signals that are 5*STD(noise) per Halouska et al. (2006)
Data=load(filename);
NoiseResonances=Data(find(Data(:,1)<-1));

for j=1:size(NoiseResonances,1)
    NoiseIntensities(j)=Absorbance_OnePoint(Data,NoiseResonances(j));
end

NoiseSTD=std(NoiseIntensities);
Noisethreshold=5*NoiseSTD;

Intensities=Data(:,2);
Intensities_corrected=Intensities;
for i=1:size(Intensities_corrected,1)
    if Intensities_corrected(i) < Noisethreshold
       Intensities_corrected(i) = 0;
    end
end
    
Data_final=[Data(:,1), Intensities_corrected];
dlmwrite([filename(1:end-4) '_denoised.txt'],Data_final,'\t');

disp(['Finished denoising ' char(filename)])
end

% Function for binning NMR spectra
function NMR_Bin(filename,binsize)

Data=load(filename);     
PPM=Data(:,1);        
Intensity=Data(:,2);  
PPM2=PPM';
Intensity2=Intensity';
PPM3=vec2mat(PPM2,binsize);
Intensity3=vec2mat(Intensity2,binsize);
PPM4=PPM3';
Intensity4=Intensity3';
PPM5=mean(PPM4);
Intensity5=mean(Intensity4);
PPM6=PPM5';
Intensity6=Intensity5';
Output=[PPM6, Intensity6];
binvalue=num2str(binsize);
dlmwrite([filename(1:end-4) '_' binvalue 'bin.txt'],Output,'\t');

disp(['Finished binning ' char(filename)])
end
 
% Function for finding intensity at specified chemical shift (same function
% for absorbance value at a wavelength)
function Abs_OnePoint = Absorbance_OnePoint(Data,Wavelength)

[~,row] = ismember(Wavelength,Data(:,1));
Abs_OnePoint = Data(row,2);

end