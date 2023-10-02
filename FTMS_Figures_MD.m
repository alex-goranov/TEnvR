function FTMS_Figures_MD (filename)
%% Description

% This code loads a list of formulas and plots a vK as well as a H/C vs MW
% plot with color-codes formulas based on mass defect (MD).

% Example:
% FTMS_Figures_MD('Sample 1_Final.xlsx')

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

% Load data and sort it by ExactMass
[Data,~]=xlsread(filename);
Data=sortrows(Data,format.Column_ExactMass); 

Data(:,format.Column_ExactMass)=Data(:,format.Column_C)*format.Mass_12C+Data(:,format.Column_H)*format.Mass_1H+Data(:,format.Column_N)*format.Mass_14N+...
    Data(:,format.Column_O)*format.Mass_16O+Data(:,format.Column_S)*format.Mass_32S+Data(:,format.Column_P)*format.Mass_31P+Data(:,format.Column_E)*format.Mass_Heteroelement;

Masses=Data(:,format.Column_ExactMass); 
NM=fix(Masses);
MD=round(Masses-NM,format.Precision);

Data_MD=[Data, NM, MD];
Data0 = Data_MD(MD > 0.95,:);
Data1 = Data_MD((0.0 <= MD) & (MD < 0.1),:);
Data2 = Data_MD((0.1 <= MD) & (MD < 0.2),:);
Data3 = Data_MD((0.2 <= MD) & (MD < 0.3),:);
Data4 = Data_MD((0.3 <= MD) & (MD < 0.4),:);
Data5 = Data_MD((0.4 <= MD) & (MD < 0.5),:);
Data6 = Data_MD((0.5 <= MD) & (MD < 0.6),:);
Data7 = Data_MD((0.6 <= MD) & (MD < 0.7),:);

% Check
if size(Data_MD,1) ~= size(Data0,1)+size(Data1,1)+size(Data2,1)+size(Data3,1)+size(Data4,1)+size(Data5,1)+size(Data6,1)+size(Data7,1)
    error('Error! Something is wrong, formula pools overlap!')
end

%% Figures

figure(1)      
    hold on
    title(filename(1:end-5),'fontsize',13,'Units','inches','Position',[2.035717129707418,3.63375,0],'Interpreter','none');
    xlabel('Nominal Mass','fontsize',13,'Units','inches','Position',[2.035716388906909,-0.24225,0]); 
    ylabel('Mass Defect','fontsize',13,'Units','inches','Position',[-0.407142873321261,1.781251698732382,0]);
    plot(Data0(:,format.Column_AImod+1),Data0(:,format.Column_AImod+2),'.','MarkerSize',8,'DisplayName','Mass Defect 0.95 - 1','Color',[0.00 0.00 0.00]); 
    plot(Data1(:,format.Column_AImod+1),Data1(:,format.Column_AImod+2),'.','MarkerSize',6,'DisplayName','Mass Defect 0 - 0.1','Color',[1.00 0.07 0.65]); 
    plot(Data2(:,format.Column_AImod+1),Data2(:,format.Column_AImod+2),'.','MarkerSize',6,'DisplayName','Mass Defect 0.1 - 0.2','Color',[0.00 0.00 1.00]); 
    plot(Data3(:,format.Column_AImod+1),Data3(:,format.Column_AImod+2),'.','MarkerSize',8,'DisplayName','Mass Defect 0.2 - 0.3','Color',[1.00 0.00 0.00]); 
    plot(Data4(:,format.Column_AImod+1),Data4(:,format.Column_AImod+2),'.','MarkerSize',11,'DisplayName','Mass Defect 0.3 - 0.4','Color',[0.10 0.74 0.10]); 
    plot(Data5(:,format.Column_AImod+1),Data5(:,format.Column_AImod+2),'.','MarkerSize',11,'DisplayName','Mass Defect 0.4 - 0.5','Color',[0.06 1.00 1.00]); 
    plot(Data6(:,format.Column_AImod+1),Data6(:,format.Column_AImod+2),'.','MarkerSize',11,'DisplayName','Mass Defect 0.5 - 0.6','Color',[0.06 1.00 1.00]); 
    plot(Data7(:,format.Column_AImod+1),Data7(:,format.Column_AImod+2),'.','MarkerSize',11,'DisplayName','Mass Defect 0.6 - 0.7','Color',[0.06 1.00 1.00]);  
    set(gca,'fontweight','bold')
    set(gca,'ycolor','k')
    set(gca,'xcolor','k')
    hold off
set(gcf,'WindowState','normal')
print(gcf,['MD_MW_' filename(1:end-5) '.png'],'-dpng','-r500');

figure(2) 
    hold on
    xlim([0 1.2]); ylim([0 2.5]);
    xlabel('O/C Ratio','fontsize',13,'Units','inches','Position',[2.035716388906909,-0.24225,0]); 
    ylabel('H/C Ratio','fontsize',13,'Units','inches','Position',[-0.407142873321261,1.781251698732382,0]);
    title(filename(1:end-5),'Color','k','fontsize',13,'Units','inches','Position',[2.035717129707418,3.63375,0],'Interpreter','none'); 
    plot(Data0(:,format.Column_OC),Data0(:,format.Column_HC),'.','MarkerSize',8,'DisplayName','Mass Defect 0.95 - 1','Color',[0.00 0.00 0.00]); 
    plot(Data1(:,format.Column_OC),Data1(:,format.Column_HC),'.','MarkerSize',6,'DisplayName','Mass Defect 0 - 0.1','Color',[1.00 0.07 0.65]); 
    plot(Data2(:,format.Column_OC),Data2(:,format.Column_HC),'.','MarkerSize',6,'DisplayName','Mass Defect 0.1 - 0.2','Color',[0.00 0.00 1.00]); 
    plot(Data3(:,format.Column_OC),Data3(:,format.Column_HC),'.','MarkerSize',8,'DisplayName','Mass Defect 0.2 - 0.3','Color',[1.00 0.00 0.00]); 
    plot(Data4(:,format.Column_OC),Data4(:,format.Column_HC),'.','MarkerSize',11,'DisplayName','Mass Defect 0.3 - 0.4','Color',[0.10 0.74 0.10]); 
    plot(Data5(:,format.Column_OC),Data5(:,format.Column_HC),'.','MarkerSize',11,'DisplayName','Mass Defect 0.4 - 0.5','Color',[0.06 1.00 1.00]); 
    plot(Data6(:,format.Column_OC),Data6(:,format.Column_HC),'.','MarkerSize',11,'DisplayName','Mass Defect 0.5 - 0.6','Color',[0.06 1.00 1.00]); 
    plot(Data7(:,format.Column_OC),Data7(:,format.Column_HC),'.','MarkerSize',11,'DisplayName','Mass Defect 0.6 - 0.7','Color',[0.06 1.00 1.00]);  
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
    hold off
set(gcf,'WindowState','normal')
print(gcf,['MD_vK_' filename(1:end-5) '.png'],'-dpng','-r500');  
disp(['Plotted Mass Defect figures on ' char(filename) '!'])
end