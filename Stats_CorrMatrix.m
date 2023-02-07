function Stats_CorrMatrix(filename_externalvariables,type)
%% Description
% This code generates a correlation matrix on a matrix of variables. Each 
% variable is correlated with each variable. Three types of correlations 
% are possible: Pearson, Kendall, and Spearman. 

% Examples: 
% Stats_CorrMatrix('ExternalVariables.xlsx','Pearson')
% Stats_CorrMatrix('ExternalVariables.xlsx','Kendall')
% Stats_CorrMatrix('ExternalVariables.xlsx','Spearman')

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

%% Stage 1: Load and organize data

clearvars -except filename_externalvariables type
close all

% External Variables Import
[~, ~, Input1] = xlsread(filename_externalvariables);
Input1_orig=Input1;
Label_ExternalVariable=Input1(1,2:end);
Input1_original=Input1;
Input1(1,:)=[];
Input1=sortrows(Input1,1); % Sort by filename

ExternalVariables=cell2mat(Input1(:,2:end));

%% Stage 2: Correlations and Data Export

if strcmp(type,'Pearson')
    [coeff,pvalues] = corr(ExternalVariables,'Type','Pearson');
    R2=coeff.^2;
    
    xlswrite('CorrMatrix_Pearson.xlsx',Input1_orig,1,'A1'); % original data = paper trail
    xlswrite('CorrMatrix_Pearson.xlsx',coeff,'Coeffs','B2');
    xlswrite('CorrMatrix_Pearson.xlsx',Label_ExternalVariable,'Coeffs','B1');
    xlswrite('CorrMatrix_Pearson.xlsx',Label_ExternalVariable','Coeffs','A2');   

    xlswrite('CorrMatrix_Pearson.xlsx',R2,'Coeffs2','B2');
    xlswrite('CorrMatrix_Pearson.xlsx',Label_ExternalVariable,'Coeffs2','B1');
    xlswrite('CorrMatrix_Pearson.xlsx',Label_ExternalVariable','Coeffs2','A2'); 

    xlswrite('CorrMatrix_Pearson.xlsx',pvalues,'pvalues','B2');
    xlswrite('CorrMatrix_Pearson.xlsx',Label_ExternalVariable,'pvalues','B1');
    xlswrite('CorrMatrix_Pearson.xlsx',Label_ExternalVariable','pvalues','A2'); 
    
    disp('Correlational analysis (Pearson) - Complete!')  
elseif strcmp(type,'Kendall')
    [coeff,pvalues] = corr(ExternalVariables,'Type','Kendall');
    R2=coeff.^2;
    
    xlswrite('CorrMatrix_Kendall.xlsx',Input1_orig,1,'A1'); % original data = paper trail
    xlswrite('CorrMatrix_Kendall.xlsx',coeff,'Coeffs','B2');
    xlswrite('CorrMatrix_Kendall.xlsx',Label_ExternalVariable,'Coeffs','B1');
    xlswrite('CorrMatrix_Kendall.xlsx',Label_ExternalVariable','Coeffs','A2');   

    xlswrite('CorrMatrix_Kendall.xlsx',R2,'Coeffs2','B2');
    xlswrite('CorrMatrix_Kendall.xlsx',Label_ExternalVariable,'Coeffs2','B1');
    xlswrite('CorrMatrix_Kendall.xlsx',Label_ExternalVariable','Coeffs2','A2'); 

    xlswrite('CorrMatrix_Kendall.xlsx',pvalues,'pvalues','B2');
    xlswrite('CorrMatrix_Kendall.xlsx',Label_ExternalVariable,'pvalues','B1');
    xlswrite('CorrMatrix_Kendall.xlsx',Label_ExternalVariable','pvalues','A2'); 
    
    disp('Correlational analysis (Kendall) - Complete!') 
elseif strcmp(type,'Spearman')
    
    for i=1:(size(ExternalVariables,2)-1)
        for j=(i+1):size(ExternalVariables,2)
            [coeff(i,j),~,pvalues(i,j)] = spear(ExternalVariables(:,i),ExternalVariables(:,j));    
        end
    end
       
    R2=coeff.^2;
    
    xlswrite('CorrMatrix_Spearman.xlsx',Input1_orig,1,'A1'); % original data = paper trail
    xlswrite('CorrMatrix_Spearman.xlsx',coeff,'Coeffs','B2');
    xlswrite('CorrMatrix_Spearman.xlsx',Label_ExternalVariable,'Coeffs','B1');
    xlswrite('CorrMatrix_Spearman.xlsx',Label_ExternalVariable','Coeffs','A2');   

    xlswrite('CorrMatrix_Spearman.xlsx',R2,'Coeffs2','B2');
    xlswrite('CorrMatrix_Spearman.xlsx',Label_ExternalVariable,'Coeffs2','B1');
    xlswrite('CorrMatrix_Spearman.xlsx',Label_ExternalVariable','Coeffs2','A2'); 

    xlswrite('CorrMatrix_Spearman.xlsx',pvalues,'pvalues','B2');
    xlswrite('CorrMatrix_Spearman.xlsx',Label_ExternalVariable,'pvalues','B1');
    xlswrite('CorrMatrix_Spearman.xlsx',Label_ExternalVariable','pvalues','A2'); 
    
    disp('Correlational analysis (Spearman) - Complete!') 
else
    error('Error! The type argument of the funciton is wrong. Choose between Pearson, Kendall, and Spearman.')
end

end