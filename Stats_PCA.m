function Stats_PCA(InputFilename,DataType)
%% Description

% This code performs principal component analysis (PCA) on a matrix of variables.

% Examples:
% Stats_PCA('UVVIS_MasterReport.xlsx','UVVIS')
% Stats_PCA('EEM_AlignmentMatrix.xlsx','EEM')
% Stats_PCA('FTMS Alignment Matrix.xlsx','FTICRMS')
% Stats_PCA('NMR Alignment Matrix.xlsx','NMR')
% Stats_PCA('ExternalVariables.xlsx','Variables')

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

%% Code

clearvars -except InputFilename DataType
close all
tic

format=FTMS_ConfigurationToolbox;

if strcmpi(DataType,'UVVIS')
    SpectralData  = true;   % data is spectral/continious;
    Variable_text = false; % variables are numeric (Wavelength)
elseif strcmpi(DataType,'EEM')
    SpectralData  = true;   % data is spectral/continious;
    Variable_text = true;  % variables are text (e.g., EX230EM400)  
elseif strcmpi(DataType,'NMR')
    SpectralData  = true;   % data is spectral/continious;
    Variable_text = false; % variables are numeric (Chemical shift)   
elseif strcmpi(DataType,'FTICRMS')
    SpectralData  = false;   % data is not spectral - discrete!
    Variable_text = false; % variables are numeric (ExactMass)  
elseif strcmpi(DataType,'Variables')
    SpectralData  = false;   % data is not spectral - discrete! 
    Variable_text = true;  % variables are text (e.g., TOC, pH, salinity,...)   
else
    error('Error!Need to set the correct DataType!')
end

OutputFilename=['PCA_Results (' DataType ').xlsx'];

%% Load matrix and obtain lists of dependent/independent variables and samples names

if strcmpi(DataType,'Variables')
    T=readtable(InputFilename,'ReadVariableNames',true,'ReadRowNames',true,'VariableNamingRule','preserve');
else
    T=readtable(InputFilename,'ReadVariableNames',true,'ReadRowNames',true,'VariableNamingRule','preserve','Sheet','Normalized');
end
Names_Variables=T.Properties.RowNames;
Matrix_Original=T.Variables;
Names_Samples=T.Properties.VariableNames;

if strcmpi(DataType,'FTICRMS')
    Names_Samples=Names_Samples(1,1:end-12);
    Matrix_Original=Matrix_Original(:,1:size(T.Properties.VariableNames,2)-12); % Cuts 12 last columns with formula info
elseif strcmpi(DataType,'Variables')
    Matrix_Original=Matrix_Original';
    Names_Samples_temp=Names_Samples;
    Names_Variables_temp=Names_Variables;
    Names_Samples=Names_Variables_temp';
    Names_Variables=Names_Samples_temp';
end

%% Check input

[Size_Variables,~]=size(Names_Variables);
[~,Size_Samples]=size(Names_Samples);
Size_Matrix=size(Matrix_Original);

if Size_Variables~= Size_Matrix(1)
    error('Error! Check the dimensions of your matrix!')
end
if Size_Samples~= Size_Matrix(2)
    error('Error! Check the dimensions of your matrix!')
end

%% Perform PCA on the Matrix

Matrix=Matrix_Original';  % Transpose the data - samples must become rows while variables become columns. 
[Matrix_Rows,Matrix_Columns]=size(Matrix);

if SpectralData % For NMR/EEM/UVVIS spectral data
    Matrix_Mean=mean(Matrix,2); 
    Matrix=Matrix-repmat(Matrix_Mean,1,Matrix_Columns); 
    Covariance=cov(Matrix);
else % For other data (FTICRMS formulas or combined variables...)
    Matrix=Matrix-ones(Matrix_Rows,1)*mean(Matrix);
    Matrix=Matrix./(ones(Matrix_Rows,1)*std(Matrix,1));
    Covariance=cov(Matrix);
end

Covariance=fillmissing(Covariance,'constant',0); 
[Loadings,Diagonal,~]=svd(Covariance);
Scores=Matrix*Loadings;

Diagonal_RowSums=sum(Diagonal,2); 
Diagonal_TotalSum=sum(Diagonal_RowSums);
Diagonal_PCVariance=(Diagonal_RowSums*100)./Diagonal_TotalSum;

if size(Diagonal_RowSums,1) > 20
    CompThreshold = 20;
else
    CompThreshold = size(Diagonal_RowSums,1);
end

% Calc p-values for noncontinuous data
if SpectralData == 0
    for i=1:CompThreshold % Calc p-values only fore the first 20 principal components
        for j=1:Size_Variables
            [~,p] = corrcoef(Scores(:,i),Matrix(:,j));
            pvalues(i,j)=p(1,2);
        end                             
    end
end

if format.pvalue_adjustment
    for i=1:size(pvalues,1) % Calc q-values for each PC set of p-values
        qvalues(i,:)=pval_adjust(pvalues(i,:),format.p_adjust_method);
    end
end

%% Export a Report
warning('off','MATLAB:xlswrite:AddSheet');

LabelDiagonal=round(Diagonal_PCVariance,0);
LabelDiagonal=string(LabelDiagonal);
LabelDiagonal=strcat(LabelDiagonal,'%');
LabelDiagonal=transpose(LabelDiagonal);

LabelPC=1:1:size(Diagonal_PCVariance,1);
LabelPC=string(LabelPC);
LabelPC=strcat('PC',' ',LabelPC,'(',LabelDiagonal,')');

% Scores
output1=[Names_Samples',string(Scores)];
output1=[' ',LabelPC;output1];

% Loadings
output2=[Names_Variables,string(Loadings)];
output2=[' ',LabelPC;output2];

% Diagonal
output3b=[transpose(LabelPC),string(Diagonal_PCVariance)];
output3a={'Component' 'Variance (%)'};
output3=[output3a;output3b];

% pvalues
if SpectralData == 0
    output4=[Names_Variables,string(pvalues')];
    if size(Diagonal_RowSums,1) < 20
        output4=[' ',LabelPC;output4]; 
    else
       output4=[' ',LabelPC(1:20);output4]; 
    end
end

% qvalues
if format.pvalue_adjustment
    output5=[Names_Variables,string(qvalues')];
    if size(Diagonal_RowSums,1) < 20
        output5=[' ',LabelPC;output5]; 
    else
       output5=[' ',LabelPC(1:20);output5]; 
    end
end

% Export everything
if size(Diagonal_RowSums,1) < 20
    if Variable_text
        xlswrite(OutputFilename,[[' ';Names_Variables],[string(char(Names_Samples))';string(Matrix_Original)]]);
    else
        xlswrite(OutputFilename,[[' ';string(Names_Variables)],[Names_Samples;string(Matrix_Original)]]);
    end
    xlswrite(OutputFilename,output1,'Scores');
    xlswrite(OutputFilename,output2,'Loadings');
    xlswrite(OutputFilename,output3,'Eigenvalues %');
    
    if SpectralData == 0
        xlswrite(OutputFilename,output4,'pvalues');
    end
    
    if format.pvalue_adjustment
        xlswrite(OutputFilename,output5,'qvalues');
    end

else % If the matrix is too big, the computer crashes during export. Export only 20 PC.
    if Variable_text
        xlswrite(OutputFilename,[[' ';Names_Variables],[string(char(Names_Samples))';string(Matrix_Original)]]);
    else
        xlswrite(OutputFilename,[[' ';string(Names_Variables)],[Names_Samples;string(Matrix_Original)]]);
    end
    xlswrite(OutputFilename,output1(:,1:21),'Scores');
    xlswrite(OutputFilename,output2(:,1:21),'Loadings');
    xlswrite(OutputFilename,output3(1:21,:),'Eigenvalues %');
    
    if SpectralData == 0
        xlswrite(OutputFilename,output4,'pvalues');
    end

    if format.pvalue_adjustment
        xlswrite(OutputFilename,output5,'qvalues');
    end

end

disp(['Principal component analysis using ' char(DataType) ' - Complete! (' num2str(toc) ' seconds)'])
end