function FTMS_KMD2(filename1,filename2,KMD)
%% Description

% This code compares two samples. It finds common formulas between the 
% two samples and finds the unique formulas of each. It then plots them on a
% KMD plot based on user-defined KMD series.

% Example:
% FTMS_KMD2('Sample 1_Final.xlsx','Sample 2_Final.xlsx','COO')

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

close all
format=FTMS_ConfigurationToolbox;

% Load data 
[Data1,~] = xlsread(filename1); 
Data2 = xlsread(filename2); 

Data1(:,format.Column_ExactMass)=Data1(:,format.Column_C)*format.Mass_12C+Data1(:,format.Column_H)*format.Mass_1H+Data1(:,format.Column_N)*format.Mass_14N+...
    Data1(:,format.Column_O)*format.Mass_16O+Data1(:,format.Column_S)*format.Mass_32S+Data1(:,format.Column_P)*format.Mass_31P+Data1(:,format.Column_E)*format.Mass_Heteroelement;

Data2(:,format.Column_ExactMass)=Data2(:,format.Column_C)*format.Mass_12C+Data2(:,format.Column_H)*format.Mass_1H+Data2(:,format.Column_N)*format.Mass_14N+...
    Data2(:,format.Column_O)*format.Mass_16O+Data2(:,format.Column_S)*format.Mass_32S+Data2(:,format.Column_P)*format.Mass_31P+Data2(:,format.Column_E)*format.Mass_Heteroelement;

Data1 = sortrows(Data1,format.Column_ExactMass); 
Data2 = sortrows(Data2,format.Column_ExactMass); 

Masses1=round(Data1(:,format.Column_ExactMass),format.Precision);
Masses2=round(Data2(:,format.Column_ExactMass),format.Precision);

% Identify Common 
Common=intersect(Masses1,Masses2);

% Find the row location of common formulas
[~,~,index1common]=intersect(Common,Masses1); 
[~,~,index2common]=intersect(Common,Masses2);

% Check
if size(index1common,1) ~= size(index2common,1)
    error('Error! You have the same molecular formula assigned to more than one peak! Check your formula lists for ExactMass duplicates!')
end

% Extract common formulas
Common1=Data1(index1common,:);
Common2=Data2(index2common,:);

% Identify Unique
index1unique=[1:1:size(Data1(:,format.Column_ExactMass),1)]'; index1unique(index1common)=[];
index2unique=[1:1:size(Data2(:,format.Column_ExactMass),1)]'; index2unique(index2common)=[];

% Extract unique formulas
Unique1=Data1(index1unique,:);
Unique2=Data2(index2unique,:);

% KMD calculations
ExactMass_PhotoLabile=Unique1(:,format.Column_ExactMass);
ExactMass_PhotoResistant=Common1(:,format.Column_ExactMass);
ExactMass_PhotoProduced=Unique2(:,format.Column_ExactMass);

% KMD calculations
if strcmp(KMD,'CH2')
     FunctionalityExactMass=format.Mass_12C+format.Mass_1H+format.Mass_1H;    
elseif strcmp(KMD,'H2')
    FunctionalityExactMass=format.Mass_1H+format.Mass_1H;
elseif strcmp(KMD,'O')
    FunctionalityExactMass=format.Mass_16O;
elseif strcmp(KMD,'CO')
    FunctionalityExactMass=format.Mass_12C+format.Mass_16O;    
elseif strcmp(KMD,'COO')
    FunctionalityExactMass=format.Mass_12C+format.Mass_16O+format.Mass_16O;
elseif strcmp(KMD,'CH2O')
     FunctionalityExactMass=format.Mass_12C+format.Mass_1H+format.Mass_1H+format.Mass_16O;
elseif strcmp(KMD,'H2O')
    FunctionalityExactMass=format.Mass_1H+format.Mass_1H+format.Mass_16O;
elseif strcmp(KMD,'NH3')
    FunctionalityExactMass=format.Mass_14N+format.Mass_1H+format.Mass_1H+format.Mass_1H;
elseif strcmp(KMD,'HE')
    FunctionalityExactMass=format.Mass_Heteroelement+format.Mass_1H;
else
    error('KMD Functionality not defined! Go to lines 17 and below and add the one you want.');
end

FunctionalityRoundMass=round(FunctionalityExactMass,0);

KM_PhotoLabile=fix(((ExactMass_PhotoLabile*FunctionalityRoundMass)/FunctionalityExactMass));
KM_PhotoResistant=fix(((ExactMass_PhotoResistant*FunctionalityRoundMass)/FunctionalityExactMass));
KM_PhotoProduced=fix(((ExactMass_PhotoProduced*FunctionalityRoundMass)/FunctionalityExactMass));
KMD_PhotoLabile=round((((ExactMass_PhotoLabile*FunctionalityRoundMass)/FunctionalityExactMass)-KM_PhotoLabile),format.Precision);
KMD_PhotoResistant=round((((ExactMass_PhotoResistant*FunctionalityRoundMass)/FunctionalityExactMass)-KM_PhotoResistant),format.Precision);
KMD_PhotoProduced=round((((ExactMass_PhotoProduced*FunctionalityRoundMass)/FunctionalityExactMass)-KM_PhotoProduced),format.Precision);

% Statistics
NumFormulas_Unique1=size(Unique1,1);
NumFormulas_Unique2=size(Unique2,1);
NumFormulas_Common=size(Common,1);
TotalFormulas=size(Data1,1)+size(Data2,1)-NumFormulas_Common;

Unique1ForPer=round(NumFormulas_Unique1*100/TotalFormulas,0);
Unique2ForPer=round(NumFormulas_Unique2*100/TotalFormulas,0);
CommonForPer=round(NumFormulas_Common*100/TotalFormulas,0);

%% Figures
figurename=[char(filename1(1:end-5)) '_' char(filename2(1:end-5))];
figure(1) % KMD series comparison
    hold on
    title(['KMD Plot (' KMD ') of ' char(filename1(1:end-5)) ' and ' char(filename2(1:end-5))],'Interpreter','none');
    xlabel(['Kendrick Nominal Mass (' KMD ')']); 
    ylabel(['Kendrick Mass Defect (' KMD ')']);
    
    String_Labile = ['Unique to Sample 1 (' num2str(NumFormulas_Unique1) ', ' num2str(Unique1ForPer) '%)' '  '];
    String_Produced = ['Unique to Sample 2 (' num2str(NumFormulas_Unique2) ', ' num2str(Unique2ForPer) '%)' '  '];
    String_Resistant = ['Common (' num2str(NumFormulas_Common) ', ' num2str(CommonForPer) '%)' '  '];       
    
    A=plot(KM_PhotoLabile,KMD_PhotoLabile,'b.','MarkerSize',13,'DisplayName',String_Labile);            %Photo/Bio-Labile RED
    B=plot(KM_PhotoProduced,KMD_PhotoProduced,'r.','MarkerSize',13,'DisplayName',String_Produced);      %Photo/Bio-Produced BLUE
    C=plot(KM_PhotoResistant,KMD_PhotoResistant,'w.','Markersize',15, 'DisplayName',String_Resistant);  %Photo/Bio-Resistant GREY
    set(C,'Color',[0.50,0.50,0.50]);
    legend([A B C]);
    h=get(gca,'Children');
    set(gca,'Children',[h(3) h(2) h(1)])
    set(gca,'fontweight','bold')
    hold off
%set(gcf,'WindowState','maximize')
print(gcf,['Compare_KMD_' KMD '_' figurename],'-dpng','-r500');

disp(['Finished comparing ' char(filename1) ' and ' char(filename2) ' using ' char(KMD) ' series.'])
end