%% Description

% This automation code can be used for:
    % 1) Automatically processing multiple FT-ICR-MS spectra (FTMS_Process)
    % 2) Automatically exploring multiple formula lists using Compound Classes (FTMS_CompoundClass)
    % 3) Automatically computing metrics for multiple formula lists (FTMS_Metrics)
    % 4) Automatically use other TEnvR functions on multiple formula lists (FTMS_Figures, FTMS_Peptides,...)

% Example:
% Run each section individually (use button "Run Section" in the Editor tab

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

%% Section 1: Automation for FTMS_Process 

%  Combined Peak Refinement, Formula Assignemnt & Formula Refinement)
clc
clear
close all

% Collect all filenames in the folder
Files=dir('*.txt'); 
Files=natsortfiles(Files);

% You will need to key in the name of the file containing blank spectrum.
% Input without apostrophes. For example: Blank PPL.txt
blankname=inputdlg('Enter the filename of the blank','Blank',1);

% Remove the blank from the structure containing the peak lists of the actual samples
Files_Filenames={Files.name};
index = find(strcmp(Files_Filenames,blankname));
Files(index)=[];

% Apply FTMS_Process on all samples
for z=1:size(Files,1)
  FileName=Files(z).name;
  FTMS_Process(FileName,char(blankname),FTMS_ConfigurationAssignment);
end

disp('Automated processing with FTMS_Process - Complete!')

%% Section 2: Automation for data exploration using FTMS_CompoundClass

clc
clear
close all

Files=dir('*_Final.xlsx'); 
Files=natsortfiles(Files);

for z=1:size(Files,1)
  FileName=Files(z).name;
  stat=FTMS_CompoundClass(FileName);
  stats(z,:)=stat;
end

Stats_TitlesColumns={'Total Peaks','CHO #','CHON #','CHO+CHON #','CHOS #','CHOP #',	'CHONS #',' CHONP #','CHOSP #','CHONSP #','Cl #','Others #','CHO num%',...
    'CHON num%','CHO+CHON num%','CHOS num%','CHOP num%','CHONS num%','CHONP num%','CHOSP num%','CHONSP num%','Cl num%','Others num%','CHO int%','CHON int%',...
    'CHO+CHON int%','CHOS int%','CHOP int%','CHONS int%','CHONP int%','CHOSP int%','CHONSP int%','Cl int%','Others int%','BC #','SCA #','Lignin #','Tannin #',...
    'Sulfonic Acid #','Amino Sugar #','Sugar #','Protein #','Fatty Acid #','Lipid #','Unsaturates #','Cl #','Extra #','BC num%','SCA num%','Lignin num%','Tannin num%',...
    'Sulfonic Acid num%','Amino Sugar num%','Sugar num%','Protein num%','Fatty Acid num%','Lipid num%','Unsaturates num%','Cl num%','Extra num%','BC int%','SCA int%',...
    'Lignin int%','Tannin int%','Sulfonic Acid int%','Amino Sugar int%','Sugar int%','Protein int%','Fatty Acid int%','Lipid int%','Unsaturates int%','Cl int%','Extra int%'};

xlswrite('FTMS CompoundClass Master Report.xlsx',Stats_TitlesColumns,1,'B1')
xlswrite('FTMS CompoundClass Master Report.xlsx',stats,1,'A2')

disp('Automated classification with FTMS_CompoundClass - Complete!')

%% Section 3: Automation for data exploration using FTMS_Metrics

clear 
close all
clc

% Import filenames
Files=dir('*_Final.xlsx');
Files=natsortfiles(Files);

for i=1:size(Files,1)
    FileName=Files(i).name;
    [AVGs,STDs,AVGsw]=FTMS_Metrics(FileName);
    Matrix_AVGs(i,:)=AVGs;
    Matrix_STDs(i,:)=STDs;
    Matrix_AVGsw(i,:)=AVGsw;
end

% Export Data
warning('off','MATLAB:xlswrite:AddSheet')
titleAVGs={'C' 'H' 'O' 'N' 'S' 'P' 'E' 'O/C' 'H/C' 'N/C' 'E/C' 'H/N' 'O/N' 'H/S' 'H/P' 'O/S' 'O/P' 'N/S' 'P/S' 'H/E' 'O/E' 'N/E' 'ExactMass' 'DBE' 'DBE/C' 'DBE/H' 'DBE/O' 'DBE-O' 'AImod' 'Ring' 'NOSC'};
titleAVGsw={'Cw' 'Hw' 'Ow' 'Nw' 'Sw' 'P' 'Ew' 'O/Cw' 'H/Cw' 'N/Cw' 'E/Cw' 'H/Nw' 'O/Nw' 'H/Sw' 'H/Pw' 'O/Sw' 'O/Pw' 'N/Sw' 'P/Sw' 'H/Ew' 'O/Ew' 'N/Ew' 'ExactMassw' 'DBEw' 'DBE/Cw' 'DBE/Hw' 'DBE/Ow' 'DBE-Ow' 'AImodw' 'Ringw' 'NOSCw'};

Files_final_Filenames={Files.name}';
xlswrite('FTMS Metrics Master Report.xlsx',Files_final_Filenames,1,'A2')
xlswrite('FTMS Metrics Master Report.xlsx',Matrix_AVGs,1,'B2')
xlswrite('FTMS Metrics Master Report.xlsx',titleAVGs,1,'B1')

xlswrite('FTMS Metrics Master Report.xlsx',Files_final_Filenames,'STDEV','A2')
xlswrite('FTMS Metrics Master Report.xlsx',Matrix_STDs,'STDEV','B2')
xlswrite('FTMS Metrics Master Report.xlsx',titleAVGs,'STDEV','B1')

xlswrite('FTMS Metrics Master Report.xlsx',Files_final_Filenames,'Metrics Weighed','A2')
xlswrite('FTMS Metrics Master Report.xlsx',Matrix_AVGsw,'Metrics Weighed','B2')
xlswrite('FTMS Metrics Master Report.xlsx',titleAVGsw,'Metrics Weighed','B1')

disp('Automated metrics computation - Complete!')

%% Section 4: Automation for generic function (FTMS_Function)

% Can be applied on FTMS_Figures, FTMS_Peptides, etc. 
clc
clear
close all

% Apply "FTMS_Function" of choice on all samples
Files=dir('*_Final.xlsx');
Files=natsortfiles(Files);

for z=1:size(Files,1)
  FileName=Files(z).name;
  FTMS_Function(FileName);
end

disp('Automated function application - Complete!')