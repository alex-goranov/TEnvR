function Stats_HCA(InputFilename,Type_Data)
%% Description

% This code performs hierarchical cluster analysis (HCA) and generates a 
% dendrogram figure. The code loads a matrix of variables. The distance 
% criterion (variable Type_Distance) and type of linkage (variable: Type_Tree)
% are adjustable:

Type_Distance='euclidean';
Type_Tree='average';

% Examples:
% Stats_HCA('UVVIS_MasterReport.xlsx','UVVIS')
% Stats_HCA('EEM_AlignmentMatrix.xlsx','EEM')
% Stats_HCA('FTMS Alignment Matrix.xlsx','FTICRMS')
% Stats_HCA('NMR Alignment Matrix.xlsx','NMR')
% Stats_HCA('ExternalVariables.xlsx','Variables')

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

%% Load matrix and obtain lists of dependent/independent variables and samples names

clearvars -except InputFilename Type_Data Type_Distance Type_Tree
close all
tic

if strcmpi(Type_Data,'Variables')
    T=readtable(InputFilename,'ReadVariableNames',true,'ReadRowNames',true,'VariableNamingRule','preserve');
else
    T=readtable(InputFilename,'ReadVariableNames',true,'ReadRowNames',true,'VariableNamingRule','preserve','Sheet','Normalized');
end

Names_Variables=T.Properties.RowNames;
Matrix_Original=T.Variables;
Names_Samples=T.Properties.VariableNames;

if strcmpi(Type_Data,'FTICRMS')
    Names_Samples=Names_Samples(1,2:end-12);
    Matrix_Original=Matrix_Original(:,1:size(Names_Samples,2));
elseif strcmpi(Type_Data,'Variables')
    Matrix_Original=Matrix_Original';
    Names_Samples_temp=Names_Samples;
    Names_Variables_temp=Names_Variables;
    Names_Samples=Names_Variables_temp;
    Names_Variables=Names_Samples_temp;
end

%% Check input

if strcmpi(Type_Data,'Variables')
    [~,Size_Variables]=size(Names_Variables);
    [Size_Samples,~]=size(Names_Samples);
    Size_Matrix=size(Matrix_Original);
else
    [Size_Variables,~]=size(Names_Variables);
    [~,Size_Samples]=size(Names_Samples);
    Size_Matrix=size(Matrix_Original);
end

if Size_Variables~= Size_Matrix(1)
    error('Error! Check the dimensions of your matrix!')
end
if Size_Samples~= Size_Matrix(2)
    error('Error! Check the dimensions of your matrix!')
end

%% Perform HCA on the Matrix

Distances=pdist(Matrix_Original',Type_Distance);
Tree=linkage(Distances,Type_Tree);

figure
    dendrogram(Tree,0);
    D=dendrogram(Tree,0,'Labels',Names_Samples,'orientation','left');
    set(D,'LineWidth',2)
    title('Dendrogram');
    xlabel('Distance'); 
    set(gca,'fontweight','bold','FontSize',12,'TickLabelInterpreter','none')
print(gcf,['HCA_' InputFilename(1:end-5) '_' Type_Data '.png'],'-dpng','-r300');

disp(['Hierarchical cluster analysis using ' char(Type_Data) ' - Complete! (' num2str(toc) ' seconds)'])
end