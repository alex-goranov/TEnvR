%% Description

% This code takes one or multiple processed EEM spectra and folds them 
% into 2D arrays. Prepares an Excel matrix for PCA or other stats.

% Example: 
% Run the code by typing EEM_Fold in the command windwow.

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
clear 

Files=dir('*_Final.csv');
Files=natsortfiles(Files);

for z=1:size(Files,1)
    FileName=Files(z).name;
    [x_fold,Int_fold]=FoldEEM(FileName);
    M(z,:)=Int_fold;
end
    
M_original=M;
x_fold_original=x_fold;
    
% replace negative and NaN values with zero
M(M<0)=0; 
M(isnan(M))=0;    
    
% Remove wavelength values with all zeros
indexes=find(sum(M == 0) == size(M,1));
M(:,indexes)=[];    
x_fold(indexes)=[];

% Normalize 
TotalSpectralInt=sum(M,2);
for i=1:size(M,1)
    RawInt=M(i,:);
    Norm(i,:)=RawInt./TotalSpectralInt(i);
end

warning('off','MATLAB:xlswrite:AddSheet')

% Export Original
xlswrite('EEM_AlignmentMatrix.xlsx',{Files.name},1,'B1')
xlswrite('EEM_AlignmentMatrix.xlsx',M_original',1,'B2')
xlswrite('EEM_AlignmentMatrix.xlsx',x_fold_original,1,'A2')

% Export Normalized
xlswrite('EEM_AlignmentMatrix.xlsx',{Files.name},'Normalized','B1')
xlswrite('EEM_AlignmentMatrix.xlsx',Norm','Normalized','B2')
xlswrite('EEM_AlignmentMatrix.xlsx',x_fold,'Normalized','A2')

disp('EEM_Folding - Complete!')

%% Internal functions

% Function for folding an individual EEM
function [x_fold,Int_fold]=FoldEEM(filename)
Data=csvread(filename);             % Import data
Ex=Data(1,:); Em=Data(:,1);         % Extract excitation and emission ranges
Ex=Ex(2:end); Em=Em(2:end);         % Extract excitation and emission data
Int=Data; Int(1,:)=[]; Int(:,1)=[]; % Extract intensity data

count=[];
    for j=1:size(Ex,2) 
        if isempty(count)
            count=0;
        end
        for z=1:size(Em,1) 
            count=count+1;
            EX_j=Ex(j);
            EM_z=Em(z);
            x_fold(count)=string(['EX' num2str(EX_j) 'EM' num2str(EM_z)]);
            Int_fold(count)=Int(z,j); % Em,Ex
        end
    end

x_fold=x_fold';
Int_fold=Int_fold';
end