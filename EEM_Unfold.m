function EEM_Unfold(filename,pcomp)
%% Description
% This code takes a report from a principal component analysis, takes  
% the loadings of the specified principal component (in 2D format) and 
% unfolds them into a regular 3D EEM format. 

% Example: 
% EEM_Unfold('PCA_Results (EEM).xlsx',1)

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

clearvars -except filename pcomp
close all

T=readtable(filename,'ReadVariableNames',true,'ReadRowNames',true,'VariableNamingRule','preserve','Sheet','Loadings');
Names_Rows=T.Properties.RowNames;
Data=T.Variables;
Loading=Data(:,pcomp);

Matrix=sep(Names_Rows,Loading);
EX=Matrix(:,1);
EM=Matrix(:,2);
Int=Matrix(:,3);

NumEX=size(unique(Matrix(:,1)),1);
NumEM=size(unique(Matrix(:,2)),1);

EX_Final = linspace(min(EX), max(EX), NumEX+1);
EM_Final = linspace(min(EM), max(EM), NumEM+1);
[X,Y] = meshgrid(EX_Final, EM_Final);
Int_Final = griddata(EX,EM,Int,X,Y);

Matrix_Final=[0,EX_Final;transpose(EM_Final),Int_Final];
csvwrite([char(filename(1:end-5)) '_comp'  num2str(pcomp) '_unfold.csv'],Matrix_Final);

disp(['Finished unfolding Principal Component ' num2str(pcomp) ' of ' char(filename)])
end


%% Internal functions

% This function separates the EX&EM values and creates a new array
function Matrix = sep(Names_Rows,Loading_2D)

for i=1:size(Names_Rows,1)
    Int(i)=Loading_2D(i);
    EX(i)=str2double(extractBetween(Names_Rows(i),"EX","EM"));
    EM(i)=str2double(extractAfter(Names_Rows(i),"EM"));
    Matrix(i,:)=[EX(i),EM(i),Int(i)];
end

end