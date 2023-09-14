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

%% Stage 1: Load and organize data

clearvars -except filename_externalvariables type
close all

global num_correl
format=FTMS_ConfigurationToolbox;

% External Variables Import
[~, ~, Input1] = xlsread(filename_externalvariables);
Input1_orig=Input1;
Label_ExternalVariable=Input1(1,2:end);
Input1_original=Input1;
Input1(1,:)=[];
Input1=sortrows(Input1,1); % Sort by filename

ExternalVariables=cell2mat(Input1(:,2:end));
num_correl=length(Label_ExternalVariable);

%% Stage 2: Correlations and Data Export

warning('off','MATLAB:xlswrite:AddSheet');

if strcmp(type,'Pearson')
    [coeff,pvalues] = corr(ExternalVariables,'Type','Pearson');
    R2=coeff.^2;

    coeff = SubWithNaN(coeff);
    pvalues = SubWithNaN(pvalues);
    R2 = SubWithNaN(R2);
    
    if format.pvalue_adjustment
        pvalues_temp = [];  % Initialize an array to store the extracted values
    
        % Extract values above the diagonal
        for i = 1:num_correl
            for j = i+1:num_correl % Start from the column after the diagonal
                pvalues_temp = [pvalues_temp, pvalues(i, j)];
            end
        end
    
        qvalues_temp=pval_adjust(pvalues_temp,format.p_adjust_method);

        % Put the extracted values back above the diagonal
        index = 1;
        for i = 1:num_correl
            for j = i+1:num_correl
                qvalues(i, j) = qvalues_temp(index);
                index = index + 1;
            end
        end
        qvalues=[qvalues;zeros(1, num_correl)];
        qvalues = SubWithNaN(qvalues);
    end
    
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
    
    if format.pvalue_adjustment
        xlswrite('CorrMatrix_Pearson.xlsx',qvalues,'qvalues','B2');
        xlswrite('CorrMatrix_Pearson.xlsx',Label_ExternalVariable,'qvalues','B1');
        xlswrite('CorrMatrix_Pearson.xlsx',Label_ExternalVariable','qvalues','A2'); 
    end
    
    disp('Correlational analysis (Pearson) - Complete!')  
elseif strcmp(type,'Kendall')
    [coeff,pvalues] = corr(ExternalVariables,'Type','Kendall');
    R2=coeff.^2;

    coeff = SubWithNaN(coeff);
    pvalues = SubWithNaN(pvalues);
    R2 = SubWithNaN(R2);
    
    if format.pvalue_adjustment
        pvalues_temp = [];  % Initialize an array to store the extracted values
    
        % Extract values above the diagonal
        for i = 1:num_correl
            for j = i+1:num_correl % Start from the column after the diagonal
                pvalues_temp = [pvalues_temp, pvalues(i, j)];
            end
        end
    
        qvalues_temp=pval_adjust(pvalues_temp,format.p_adjust_method);

        % Put the extracted values back above the diagonal
        index = 1;
        for i = 1:num_correl
            for j = i+1:num_correl
                qvalues(i, j) = qvalues_temp(index);
                index = index + 1;
            end
        end
        qvalues=[qvalues;zeros(1, num_correl)];
        qvalues = SubWithNaN(qvalues);
    end

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
    
    if format.pvalue_adjustment
        xlswrite('CorrMatrix_Kendall.xlsx',qvalues,'qvalues','B2');
        xlswrite('CorrMatrix_Kendall.xlsx',Label_ExternalVariable,'qvalues','B1');
        xlswrite('CorrMatrix_Kendall.xlsx',Label_ExternalVariable','qvalues','A2'); 
    end

    disp('Correlational analysis (Kendall) - Complete!') 

elseif strcmp(type,'Spearman')
    
    for i=1:(size(ExternalVariables,2)-1)
        for j=(i+1):size(ExternalVariables,2)
            [coeff(i,j),~,pvalues(i,j)] = spear(ExternalVariables(:,i),ExternalVariables(:,j));    
        end
    end
    
    coeff=[coeff;zeros(1, num_correl)];
    pvalues=[pvalues;zeros(1, num_correl)];
    R2=coeff.^2;

    coeff = SubWithNaN(coeff);
    pvalues = SubWithNaN(pvalues);
    R2 = SubWithNaN(R2);
    
    if format.pvalue_adjustment
        pvalues_temp = [];  % Initialize an array to store the extracted values
    
        % Extract values above the diagonal
        for i = 1:num_correl
            for j = i+1:num_correl % Start from the column after the diagonal
                pvalues_temp = [pvalues_temp, pvalues(i, j)];
            end
        end
    
        qvalues_temp=pval_adjust(pvalues_temp,format.p_adjust_method);

        % Put the extracted values back above the diagonal
        index = 1;
        for i = 1:num_correl
            for j = i+1:num_correl
                qvalues(i, j) = qvalues_temp(index);
                index = index + 1;
            end
        end
        qvalues=[qvalues;zeros(1, num_correl)];
        qvalues = SubWithNaN(qvalues);
    end
    
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

    if format.pvalue_adjustment
        xlswrite('CorrMatrix_Spearman.xlsx',qvalues,'qvalues','B2');
        xlswrite('CorrMatrix_Spearman.xlsx',Label_ExternalVariable,'qvalues','B1');
        xlswrite('CorrMatrix_Spearman.xlsx',Label_ExternalVariable','qvalues','A2'); 
    end

    disp('Correlational analysis (Spearman) - Complete!') 
else
    error('Error! The type argument of the funciton is wrong. Choose between Pearson, Kendall, and Spearman.')
end

end

%% Internal Function - substitue diagonal and below with NaN
function A  = SubWithNaN(A)

global num_correl

    A(logical(eye(size(A)))) = NaN;
    for i = 1:num_correl
        for j = 1:num_correl
            if i > j
                A(i, j) = NaN;
            end
        end
    end
end