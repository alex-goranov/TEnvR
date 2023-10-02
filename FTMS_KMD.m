function FTMS_KMD(filename,KMD_type)
%% Description

% This code plots KMD figures using user-defined KMD series. 

% Example:
% FTMS_KMD('Sample 1_Final.xlsx','COO')

%% Copyright and License Notices: 

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

% Load data and sort it by ExactMass
[Data,TXT,~]=xlsread(filename);
Headers=TXT(1,:);
Data=sortrows(Data,format.Column_ExactMass); 

Data(:,format.Column_ExactMass)=Data(:,format.Column_C)*format.Mass_12C+Data(:,format.Column_H)*format.Mass_1H+Data(:,format.Column_N)*format.Mass_14N+...
    Data(:,format.Column_O)*format.Mass_16O+Data(:,format.Column_S)*format.Mass_32S+Data(:,format.Column_P)*format.Mass_31P+Data(:,format.Column_E)*format.Mass_Heteroelement;

% KMD calculations
if strcmp(KMD_type,'CH2')
     FunctionalityExactMass=format.Mass_12C+format.Mass_1H+format.Mass_1H;    
elseif strcmp(KMD_type,'H2')
    FunctionalityExactMass=format.Mass_1H+format.Mass_1H;
elseif strcmp(KMD_type,'O')
    FunctionalityExactMass=format.Mass_16O;
elseif strcmp(KMD_type,'CO')
    FunctionalityExactMass=format.Mass_12C+format.Mass_16O;    
elseif strcmp(KMD_type,'COO')
    FunctionalityExactMass=format.Mass_12C+format.Mass_16O+format.Mass_16O;
elseif strcmp(KMD_type,'CH2O')
     FunctionalityExactMass=format.Mass_12C+format.Mass_1H+format.Mass_1H+format.Mass_16O;
elseif strcmp(KMD_type,'H2O')
    FunctionalityExactMass=format.Mass_1H+format.Mass_1H+format.Mass_16O;
elseif strcmp(KMD_type,'NH3')
    FunctionalityExactMass=format.Mass_14N+format.Mass_1H+format.Mass_1H+format.Mass_1H;
elseif strcmp(KMD_type,'NH')
    FunctionalityExactMass=format.Mass_14N+format.Mass_1H;
elseif strcmp(KMD_type,'SO3')
    FunctionalityExactMass=format.Mass_32S+format.Mass_16O+format.Mass_16O+format.Mass_16O;
elseif strcmp(KMD_type,'SO2')
    FunctionalityExactMass=format.Mass_32S+format.Mass_16O+format.Mass_16O;
elseif strcmp(KMD_type,'HE')
    FunctionalityExactMass=format.Mass_Heteroelement+format.Mass_1H;
else
    error('KMD Functionality not defined! Go to lines 17 and below and add the one you want.');
end

Masses=Data(:,format.Column_ExactMass); 
FunctionalityRoundMass=round(FunctionalityExactMass,0);
KM=fix(((Masses*FunctionalityRoundMass)/FunctionalityExactMass));
KMD=round((((Masses*FunctionalityRoundMass)/FunctionalityExactMass)-KM),format.Precision);

Data=[Data, KM, KMD];

% Identify formulas not sitting on the series ("Unique"). 
% Part of series = "duplicates" (from having duplicate or more "same" KMD values)
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

Index_Products=[1:1:size(Data_Duplicates(:,size(Data_Duplicates,2)),1)]'; Index_Products(Index_apo)=[];
Data_Products=Data_Duplicates(Index_Products,:);

% Statistics
NumFormulas_Unique=size(Data_Unique,1);         % Not on series
%NumFormulas_OnSeries=size(Data_duplicates,1);  % Total number of formulas on series
NumFormulas_Substrate=size(Data_Substrates,1);
NumFormulas_Products=size(Data_Products,1);

Data_Unique=sortrows(Data_Unique,format.Column_ExactMass); 
Data_Duplicates=sortrows(Data_Duplicates,format.Column_ExactMass); 
Data_Substrates=sortrows(Data_Substrates,format.Column_ExactMass); 
Data_Products=sortrows(Data_Products,format.Column_ExactMass); 

KM_Unique=Data_Unique(:,size(Data_Unique,2)-1);
KMD_Unique=Data_Unique(:,size(Data_Unique,2));

KM_Substrates=Data_Substrates(:,size(Data_Substrates,2)-1);
KMD_Substrates=Data_Substrates(:,size(Data_Substrates,2));

KM_Products=Data_Products(:,size(Data_Products,2)-1);
KMD_Products=Data_Products(:,size(Data_Products,2));

PerUnique=round(NumFormulas_Unique*100/size(Data,1),0);
PerSubstrate=round(NumFormulas_Substrate*100/size(Data,1),0);
PerProducts=round(NumFormulas_Products*100/size(Data,1),0);

if PerUnique+PerSubstrate+PerProducts < 98
    error('Error! The percentages do not equal 100! Revisit the calculations!')
end

%% Export
    warning('off','MATLAB:xlswrite:AddSheet')

    xlswrite(['KMD_' KMD_type '_' char(filename(1:end-5)) '.xlsx'],Headers,'Sheet1','A1')
    xlswrite(['KMD_' KMD_type '_' char(filename(1:end-5)) '.xlsx'],Data,'Sheet1','A2')

    xlswrite(['KMD_' KMD_type '_' char(filename(1:end-5)) '.xlsx'],[Headers 'KNM' 'KMD'],'Unique','A1')
    xlswrite(['KMD_' KMD_type '_' char(filename(1:end-5)) '.xlsx'],Data_Unique,'Unique','A2')

    xlswrite(['KMD_' KMD_type '_' char(filename(1:end-5)) '.xlsx'],[Headers 'KNM' 'KMD'],'Apo','A1')
    xlswrite(['KMD_' KMD_type '_' char(filename(1:end-5)) '.xlsx'],Data_Substrates,'Apo','A2')

    xlswrite(['KMD_' KMD_type '_' char(filename(1:end-5)) '.xlsx'],[Headers 'KNM' 'KMD'],'Sequential','A1')
    xlswrite(['KMD_' KMD_type '_' char(filename(1:end-5)) '.xlsx'],Data_Products,'Sequential','A2')

%% Figures

figurename=[char(filename(1:end-5)) '.png'];
figure(1)      
    hold on
    String_Unique = ['Not on KMD series (' num2str(NumFormulas_Unique) ', ' num2str(PerUnique) '%)' '  '];
    String_Substrate = ['Apo (' num2str(NumFormulas_Substrate) ', ' num2str(PerSubstrate) '%)' '  '];
    String_Products = ['Sequential (' num2str(NumFormulas_Products) ', ' num2str(PerProducts) '%)' '  '];
    title(filename(1:end-5),'fontsize',13,'Units','inches','Position',[2.035717129707418,3.63375,0],'Interpreter','none');
    xlabel(['Kendrick Nominal Mass (' KMD_type ')'],'fontsize',13,'Units','inches','Position',[2.035716388906909,-0.24225,0]); 
    ylabel(['Kendrick Mass Defect (' KMD_type ')'],'fontsize',13,'Units','inches','Position',[-0.407142873321261,1.781251698732382,0]);
    B=plot(KM_Unique,KMD_Unique,'r.','MarkerSize',11,'DisplayName',String_Unique); 
    C=plot(KM_Products,KMD_Products,'b.','MarkerSize',12,'DisplayName',String_Products); 
    A=plot(KM_Substrates,KMD_Substrates,'b.','MarkerSize',12,'DisplayName',String_Substrate); 
    set(B,'Color',[0.50,0.50,0.50]);
    set(A,'Color',[0.00,0.27,0.00]);
    set(C,'Color',[0.1,0.74,0.1]);
    H=legend([B A C]);
    set(gca,'fontweight','bold')
    set(gca,'ycolor','k')
    set(gca,'xcolor','k')
    set(gca,'OuterPosition',[0.012520737327189,0,0.9036866359447,1]);
    set(gca,'Position',[0.13,0.11,0.700357142857143,0.815]);
    set(gca,'PlotBoxAspectRatio',[1,0.874999965230625,0.874999965230625]);
    H.Units = 'inches';
    H.Position = [3.3,0.614062502483526,2.333333279627066,0.552083318432172];
    H.Orientation = 'horizontal';
    H.NumColumns = 1;
    hold off
    set(gcf,'WindowState','normal')
    print(gcf,['KMD_' KMD_type '_' figurename],'-dpng','-r500');
figure(2) 
    hold on
    xlim([0 1.2]); ylim([0 2.5]);
    xlabel('O/C Ratio','fontsize',13,'Units','inches','Position',[2.035716388906909,-0.24225,0]); 
    ylabel('H/C Ratio','fontsize',13,'Units','inches','Position',[-0.407142873321261,1.781251698732382,0]);
    title(filename(1:end-5),'Color','k','fontsize',13,'Units','inches','Position',[2.035717129707418,3.63375,0],'Interpreter','none'); 
    B=plot(Data_Unique(:,format.Column_OC),Data_Unique(:,format.Column_HC),'r.','MarkerSize',11,'DisplayName',String_Unique); 
    C=plot(Data_Products(:,format.Column_OC),Data_Products(:,format.Column_HC),'b.','MarkerSize',12,'DisplayName',String_Products);
    A=plot(Data_Substrates(:,format.Column_OC),Data_Substrates(:,format.Column_HC),'b.','MarkerSize',12,'DisplayName',String_Substrate); 
    set(B,'Color',[0.50,0.50,0.50]);
    set(A,'Color',[0.00,0.27,0.00]);
    set(C,'Color',[0.1,0.74,0.1]);
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
    print(gcf,['KMDvK_' KMD_type '_' figurename],'-dpng','-r500');  
disp(['Performed KMD analysis on ' char(filename) ' using ' char(KMD_type) ' series.'])
end