function FTMS_Peptides(filename,minAAnumber,maxAAnumber)
%% Description 

% This code loads a list of formulas and evaluates if any of the formulas are of oligopeptides

% Example
% FTMS_Peptides('Sample1_Final.xlsx',2,5)

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

clearvars -except filename minAAnumber maxAAnumber
format=FTMS_ConfigurationToolbox;

[Data,TXT] = xlsread(filename);

Data(:,format.Column_ExactMass)=Data(:,format.Column_C)*format.Mass_12C+Data(:,format.Column_H)*format.Mass_1H+Data(:,format.Column_N)*format.Mass_14N+...
    Data(:,format.Column_O)*format.Mass_16O+Data(:,format.Column_S)*format.Mass_32S+Data(:,format.Column_P)*format.Mass_31P+Data(:,format.Column_E)*format.Mass_Heteroelement;
Data(:,format.Column_ExactMass)=round(Data(:,format.Column_ExactMass),format.Precision+1); % Round to the uncertainty in the original data

Data=sortrows(Data,format.Column_ExactMass);
FormulaType = string(TXT(2:end,format.Column_Type));

AA.mass(1)  = 89.047678473;  AA.name(1)  = 'A';
AA.mass(2)  = 174.111675712; AA.name(2)  = 'R';
AA.mass(3)  = 132.053492132; AA.name(3)  = 'N';
AA.mass(4)  = 133.037507717; AA.name(4)  = 'D';
AA.mass(5)  = 147.053157781; AA.name(5)  = 'E'; 
AA.mass(6)  = 146.069142196; AA.name(6)  = 'Q';
AA.mass(7)  = 75.032028409;  AA.name(7)  = 'G';
AA.mass(8)  = 155.069476547; AA.name(8)  = 'H';
AA.mass(9)  = 131.058243159; AA.name(9)  = 'O';
AA.mass(10) = 131.094628665; AA.name(10) = 'L'; % Or I, they are structural isomers!
AA.mass(11) = 146.105527702; AA.name(11) = 'K';
AA.mass(12) = 165.078978601; AA.name(12) = 'F';
AA.mass(13) = 115.063328537; AA.name(13) = 'P';
AA.mass(14) = 129.042593095; AA.name(14) = 'U';
AA.mass(15) = 105.042593095; AA.name(15) = 'S';
AA.mass(16) = 119.058243159; AA.name(16) = 'T';
AA.mass(17) = 204.089877638; AA.name(17) = 'W';
AA.mass(18) = 181.073893223; AA.name(18) = 'Y';
AA.mass(19) = 117.078978601; AA.name(19) = 'V';
AA.mass(20) = 132.089877638; AA.name(20) = 'X'; % Ornitine, unofficial label!

i = minAAnumber:maxAAnumber;
W = {};
p = 0;
for j = i 
    p = p+1;
   
    r = combn([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20],j);
    
    kk = unique(sort(r,2), 'rows');
    
    W{p,1} = kk;
end

mH2O=1.007825032 + 1.007825032 + 15.994914622;
for k = 1 : size(W,1)
    
    a = AA.name(W{k});
    b = AA.mass(W{k});
    
    PeptideNames{k,1} = a;
    PeptideMasses{k,1} = b;
    
    PeptideMassesSum{k,1} = sum(PeptideMasses{k},2);
    PeptideMassesSum{k,1} = PeptideMassesSum{k,1}-k*mH2O; % H2O is removed when peptidic bond is formed
end

for z = 1:size(PeptideNames,1)  
           
    l = ' ';
    l = repmat(l, size(PeptideNames{z,1},1),(maxAAnumber-size(PeptideNames{z,1},2)));
    PeptideNames{z,1} = [PeptideNames{z,1} l];
        
end

PeptideNamesMat = cell2mat(PeptideNames);
PeptideNamesMatMasses = cell2mat(PeptideMassesSum);
PeptideNamesMatMasses = round(PeptideNamesMatMasses,format.Precision);% Round to the certainty (eliminate the uncertainty!)

% Load exact mass
Masses = Data(:,format.Column_ExactMass);
Masses = round(Masses,format.Precision); % Round to the certainty (eliminate the uncertainty!)

for q = 1:size(Masses,1)
    diff = Masses(q)-PeptideNamesMatMasses;
    error =(diff./PeptideNamesMatMasses)*10^6;
    
    % Identify identical m/z:
    location = find(error >= -format.MassAccuracy & error <= format.MassAccuracy);
    
    if ~isempty(location) 
        
        A = repmat(Data(q,:),size(location,1),1);           % A is original data
        B = repmat(size(location,1),size(location,1),1);    % B is number of possibilities
        C = PeptideNamesMat(location,:);                    % D identifies correct formulas
        D = error(location);                                % F is error on each formula
    else
        A = repmat(Data(q,:),1,1);                          % A is original data
        B = 0;                                              % B is number of possibilities
        C = repmat(' ', 1, maxAAnumber);                    % D identifies correct formulas
        D = 1000;                                           % F is error on each formula
    end
    found{q,1} = [A];
    found{q,2} = [B];
    found{q,3} = [C];
    found{q,4} = [D];
end

file1=cell2mat(found(:,1));
file2=cell2mat(found(:,2));
file3=cell2mat(found(:,3));
file4=cell2mat(found(:,4));

labels1 = TXT(1,:);
labels2 = [labels1 'Options' 'Sequence' 'Error (ppm)'];
output=[string(file1),string(file2),string(file3),string(file4)];

xlswrite([filename(1:end-5) '_Peptides.xlsx'], labels2, 1, 'A1');
xlswrite([filename(1:end-5) '_Peptides.xlsx'], output, 1, 'A2');
xlswrite([filename(1:end-5) '_Peptides.xlsx'], FormulaType, 1, [char(format.Column_Type+'A'-1) '2']);

disp(['Finished finding oligopeptides in ' char(filename)])
end