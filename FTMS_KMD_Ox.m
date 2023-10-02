function FTMS_KMD_Ox(filename)
%% Description

% %This code determines oxygenation substrates and products from 
% O, CO, and COO KMD series!

% Example:
% FTMS_KMD_Ox('Sample 1_Final.xlsx')

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
format=FTMS_ConfigurationToolbox;
global Precision
Precision=format.Precision;

% Load data and sort it by ExactMass
[Data,TXT,~]=xlsread(filename);
Headers=TXT(1,:);
Data=sortrows(Data,format.Column_ExactMass); 

Data(:,format.Column_ExactMass)=Data(:,format.Column_C)*format.Mass_12C+Data(:,format.Column_H)*format.Mass_1H+Data(:,format.Column_N)*format.Mass_14N+...
    Data(:,format.Column_O)*format.Mass_16O+Data(:,format.Column_S)*format.Mass_32S+Data(:,format.Column_P)*format.Mass_31P+Data(:,format.Column_E)*format.Mass_Heteroelement;

%% KMD and Oxygenation analysis

Masses=Data(:,format.Column_ExactMass); 

[~,index_dupl] = unique(Data(:,format.Column_ExactMass),'stable');
if index_dupl ~=size(Data,1)
    error('Error! You have multiple formulas assigned per one peak!')
end

FunctionalityExactMass_O=format.Mass_16O;
FunctionalityExactMass_COO=format.Mass_12C+format.Mass_16O+format.Mass_16O;
FunctionalityExactMass_CO=format.Mass_12C+format.Mass_16O;

[~,Data_SubstratesO,Data_productsO] = KMDanalysis(Data,Masses,FunctionalityExactMass_O);
[~,Data_SubstratesCO,Data_productsCO] = KMDanalysis(Data,Masses,FunctionalityExactMass_CO);
[~,Data_SubstratesCOO,Data_productsCOO] = KMDanalysis(Data,Masses,FunctionalityExactMass_COO);

% Combine all
Substrates=[Data_SubstratesO;Data_SubstratesCO;Data_SubstratesCOO];
Products=[Data_productsO;Data_productsCO;Data_productsCOO];

Products=Products(:,1:(size(Products,2)-2));             % Trim KM and KMD columns
Substrates=Substrates(:,1:(size(Substrates,2)-2));       % Trim KM and KMD columns

Products=sortrows(Products,format.Column_ExactMass);     % Sort based on ExactMass
Substrates=sortrows(Substrates,format.Column_ExactMass); % Sort based on ExactMass

% Removal of multiples (duplicate formulas)
[~,index_dupl] = unique(Substrates(:,format.Column_ExactMass),'stable');
Substrates = Substrates(index_dupl,:);

[~,index_dupl] = unique(Products(:,format.Column_ExactMass),'stable');
Products = Products(index_dupl,:);

%% Separate Coommon and Unique
Data(:,format.Column_ExactMass)=round(Data(:,format.Column_ExactMass),Precision);
Substrates(:,format.Column_ExactMass)=round(Substrates(:,format.Column_ExactMass),Precision);
Products(:,format.Column_ExactMass)=round(Products(:,format.Column_ExactMass),Precision);

% Common between data and substrates - remove common from Data

Common=intersect(Data(:,format.Column_ExactMass),Substrates(:,format.Column_ExactMass));
[~,~,index1c]=intersect(Common,Data(:,format.Column_ExactMass)); % index is the location of the common formulas!

index1u=[1:1:length(Data(:,format.Column_ExactMass))]'; index1u(index1c)=[];
Data1_trim=Data(index1u,:);

% Common between data and products - remove common from Data

Common=intersect(Data1_trim(:,format.Column_ExactMass),Products(:,format.Column_ExactMass));
[~,~,index1c]=intersect(Common,Data1_trim(:,format.Column_ExactMass)); % index is the location of the common formulas!

index1u=[1:1:length(Data1_trim(:,format.Column_ExactMass))]'; index1u(index1c)=[];
Unique=Data1_trim(index1u,:);

% Common between substrates and products - remove common from Data

Common=intersect(Substrates(:,format.Column_ExactMass),Products(:,format.Column_ExactMass));
[~,~,index1c]=intersect(Common,Substrates(:,format.Column_ExactMass)); % index is the location of the common formulas!

index1u=[1:1:length(Substrates(:,format.Column_ExactMass))]'; index1u(index1c)=[];
Substrates=Substrates(index1u,:);

%% Statistics

NumFormulas_Unique=size(Unique,1); %Not on series
NumFormulas_Substrate=size(Substrates,1);
NumFormulas_Products=size(Products,1);

if NumFormulas_Unique+NumFormulas_Substrate+NumFormulas_Products < size(Data,1)
    error('Error! The percentages do not equal 100! Revisit the calculations!')
end

PerUnique=round(NumFormulas_Unique*100/(NumFormulas_Unique+NumFormulas_Substrate+NumFormulas_Products));
PerSubstrate=round(NumFormulas_Substrate*100/(NumFormulas_Unique+NumFormulas_Substrate+NumFormulas_Products));
PerProducts=round(NumFormulas_Products*100/(NumFormulas_Unique+NumFormulas_Substrate+NumFormulas_Products));

if PerUnique+PerSubstrate+PerProducts < 98
    error('Error! The percentages do not equal 100! Revisit the calculations!')
end

%% Export
    warning('off','MATLAB:xlswrite:AddSheet')

    xlswrite(['KMD_Oxidation_'  char(filename(1:end-5)) '.xlsx'],Headers,'Sheet1','A1')
    xlswrite(['KMD_Oxidation_'  char(filename(1:end-5)) '.xlsx'],Data,'Sheet1','A2')

    xlswrite(['KMD_Oxidation_'  char(filename(1:end-5)) '.xlsx'],Headers,'Unique','A1')
    xlswrite(['KMD_Oxidation_'  char(filename(1:end-5)) '.xlsx'],Unique,'Unique','A2')

    xlswrite(['KMD_Oxidation_'  char(filename(1:end-5)) '.xlsx'],Headers,'Substrates','A1')
    xlswrite(['KMD_Oxidation_'  char(filename(1:end-5)) '.xlsx'],Substrates,'Substrates','A2')

    xlswrite(['KMD_Oxidation_'  char(filename(1:end-5)) '.xlsx'],Headers,'Ox products','A1')
    xlswrite(['KMD_Oxidation_'  char(filename(1:end-5)) '.xlsx'],Products,'Ox products','A2')

%% Figures
figurename=[char(filename(1:end-5)) '.png'];
figure(1) 
        hold on
        xlim([0 1.2]); ylim([0 2.5]);
        xlabel('O/C Ratio','fontsize',13,'Units','inches','Position',[2.035716388906909,-0.24225,0]); 
        ylabel('H/C Ratio','fontsize',13,'Units','inches','Position',[-0.407142873321261,1.781251698732382,0]);
        title(filename(1:end-5),'Color','k','fontsize',13,'Units','inches','Position',[2.035717129707418,3.63375,0],'Interpreter','none'); 
        String_Unique = ['Not on KMD series (' num2str(NumFormulas_Unique) ', ' num2str(PerUnique) '%)' '  '];
        String_Substrate = ['Substrates (' num2str(NumFormulas_Substrate) ', ' num2str(PerSubstrate) '%)' '  '];
        String_Products = ['Oxygenation products (' num2str(NumFormulas_Products) ', ' num2str(PerProducts) '%)' '  '];
        B=plot(Unique(:,format.Column_OC),Unique(:,format.Column_HC),'r.','MarkerSize',11,'DisplayName',String_Unique); 
        C=plot(Products(:,format.Column_OC),Products(:,format.Column_HC),'b.','MarkerSize',12,'DisplayName',String_Products);
        A=plot(Substrates(:,format.Column_OC),Substrates(:,format.Column_HC),'b.','MarkerSize',12,'DisplayName',String_Substrate); 
        set(B,'Color',[0.50,0.50,0.50]);
        set(A,'Color',[0.00,0.27,0.00]);
        set(C,'Color',[0.00,0.61,0.00]);
        H=legend([B A C]);
        h=get(gca,'Children');
        h1=refline(-1,2);           % AImod=0
        h2=refline(-0.3845,1.0584); % AImod=0.5
        h3=refline(-0.3740,0.7551); % AImod=0.67
        h1.Color='k';h2.Color='k';h3.Color='k';
        h1.LineWidth=1.5;h2.LineWidth=1.5;h3.LineWidth=1.5;
        h1.HandleVisibility='off';h2.HandleVisibility='off';h3.HandleVisibility='off';
        text(1.205, 0.77, 'AI_{MOD} = 0.00','fontweight','bold','fontsize',11)
        text(1.205, 0.57, 'AI_{MOD} = 0.50','fontweight','bold','fontsize',11)
        text(1.205, 0.28, 'AI_{MOD} = 0.67','fontweight','bold','fontsize',11)
        set(gca,'fontweight','bold')
        set(gca,'ycolor','k')
        set(gca,'xcolor','k')
        set(gca,'OuterPosition',[0.012520737327189,0,0.9036866359447,1]);
        set(gca,'Position',[0.13,0.11,0.700357142857143,0.815]);
        set(gca,'PlotBoxAspectRatio',[1,0.874999965230625,0.874999965230625]);
        H.Units = 'inches';    
        H.Position = [3.55,3.65,1.88541662630936,0.380208323399226];
        H.Orientation = 'horizontal';
        H.NumColumns = 1;
        hold off
        set(gcf,'WindowState','normal')
        print(gcf,['KMDvK_Oxidation_' figurename],'-dpng','-r500');

disp(['Performed KMD oxidation analysis on ' char(filename) ' using O, CO, and COO series.'])
end

%% Internal Functions

% Function for identifying formulas not sitting on the series ("Unique") as well as 
% formulas part of series = "duplicates" (from having duplicate or more same KMD values)
function [Data_Unique,Data_Substrates,Data_Products] = KMDanalysis(Data,Mass,FunctionalityExactMass)

global Precision

FunctionalityRoundMass=round(FunctionalityExactMass,0);
KM=fix(((Mass*FunctionalityRoundMass)/FunctionalityExactMass));
KMD=round((((Mass*FunctionalityRoundMass)/FunctionalityExactMass)-KM),Precision);

Data=[Data, KM, KMD];

[n, bin] = histc(KMD, unique(KMD));
multiple = find(n > 1);
Index_Duplicates = find(ismember(bin, multiple));
Index_Unique     = find(~ismember(bin, multiple));

Data_Duplicates=Data(Index_Duplicates,:);
Data_Unique=Data(Index_Unique,:);

% Check if the difference in KM equals the FunctionalityRoundMass % (e.g., for COO, delta KM must be = 44)

Data_Duplicates=sortrows(Data_Duplicates,size(Data_Duplicates,2));
UniqueKMDvalues=unique(Data_Duplicates(:,size(Data_Duplicates,2))); % identify unique KMD values

Data_Duplicates_verified=[];
Data_Unique_new=[];
for i=1:size(UniqueKMDvalues,1) % for each unique KMD value
    KMDvalue=UniqueKMDvalues(i);
    index=find(Data_Duplicates(:,size(Data_Duplicates,2)) == KMDvalue); 
    FormulasWithThisKMD=Data_Duplicates(index,:);
    dx = diff(FormulasWithThisKMD(:,size(Data_Duplicates,2)-1));
    
    if  all(mod(dx/FunctionalityRoundMass,1) == 0)  % all formulas are with the same correct delta KM
        Data_Duplicates_verified=[Data_Duplicates_verified;FormulasWithThisKMD];  
    else                                            % some formulas have a bad delta KM
        Data_Unique_new=[Data_Unique_new;FormulasWithThisKMD];
    end
end
    
Data_Duplicates=Data_Duplicates_verified;
Data_Unique=[Data_Unique_new;Data_Unique];

% Identify first and following members of KMD series (substrate + oxidation products) 
Data_Duplicates=sortrows(Data_Duplicates,size(Data_Duplicates,2));
[~,Index_apo,~]=unique(Data_Duplicates(:,size(Data_Duplicates,2)));
Data_Substrates=Data_Duplicates(Index_apo,:);

Index_Products=[1:1:length(Data_Duplicates(:,size(Data_Duplicates,2)))]'; Index_Products(Index_apo)=[];
Data_Products=Data_Duplicates(Index_Products,:);
end