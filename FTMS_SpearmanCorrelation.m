function FTMS_SpearmanCorrelation(filename_FTMSdata,filename_externalvariables)
%% Description

% This code performs Spearman correlations on individual formulas from MS
% with external parameters. 

% Example
% FTMS_SpearmanCorrelation('FTMS Alignment Matrix.xlsx','ExternalVariables.xlsx')

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

%% Configuration

clearvars -except filename_FTMSdata filename_externalvariables
close all
tic

filename_results='FTMS_SpearmanResults.xlsx';
format=FTMS_ConfigurationToolbox;

%% Stage 1: Load and organize external variables and FTMS data

% External Variables Import
[~, ~, Input1] = xlsread(filename_externalvariables);
Label_ExternalVariable=Input1(1,2:end);
Input1_original=Input1;
Input1(1,:)=[];
Input1=sortrows(Input1,1); % Sort by filename

Label_ExternalSamples = Input1(:,1);
ExternalVariables=cell2mat(Input1(:,2:end));

% FTMS Data import
[~,~,Input2]=xlsread(filename_FTMSdata,'Normalized Common');
Label_FTMS=Input2(1,:);

Label_FTMSSamples=Input2(1,2:end-12);
FTMS_Formulas=Input2(:,size(Label_FTMSSamples,2)+2:end);
FTMS_Formulas_Labels=[Input2(1,1),FTMS_Formulas(1,:)];
FTMS_Formulas_Data=cell2mat(FTMS_Formulas(2:end,:));
FTMS_Formulas_Data=[cell2mat(Input2(2:end,1)),FTMS_Formulas_Data];

FTMSData = Input2(:,2:size(Label_FTMSSamples,2)+1);
FTMSData=FTMSData';
FTMSData=sortrows(FTMSData,1); % Sort by filename
FTMSData = FTMSData';
Label_FTMSSamples=FTMSData(1,:);

% Check number of samples
if size(Label_ExternalSamples,1) ~= size(Label_FTMSSamples,2)
    error('Error! Not equal number of samples in the Alignment Matrix and in External Variables! Check input!')
end

% Check data alignment
for i=1:size(Label_FTMSSamples,2)
    if ~strcmpi(Label_ExternalSamples(i)',Label_FTMSSamples(i))
        error('Error! Sample names in Alignment matrix and External Variables do not match! Check input!')
    end
end

FTMS_Int=cell2mat(FTMSData(2:end,:));

% Export original data in Excel
warning('off','MATLAB:xlswrite:AddSheet')
xlswrite(filename_results,Input2,1);
xlswrite(filename_results,Input1_original,'ExternalVar');

%% Stage 2: Spearman Correlations between the FTMS data and all external parameters 

for ExternalVariable_Column=1:size(Label_ExternalVariable,2)

ExternalVariable_Label=char(Label_ExternalVariable(ExternalVariable_Column));

for i=1:size(FTMS_Int,1)
X = FTMS_Int(i,:);
X=X';
Y = ExternalVariables(:,ExternalVariable_Column);
   [r(i),t(i),p(i)] = spear(X,Y); 
end

if format.pvalue_adjustment
    q=pval_adjust(p,format.p_adjust_method);
    q(q >= 1-format.CL_alpha/100) = NaN;
    SpearmanResults=[r',t',p',q'];
    FTMS_CorrelationResults=[FTMS_Formulas_Data, SpearmanResults];
    FTMS_CorrelationResults(any(isnan(FTMS_CorrelationResults), 2), :) = [];
    
    if ~isempty(FTMS_CorrelationResults)
        % Export data in the Excel file
        xlswrite(filename_results,[FTMS_Formulas_Labels 'r-value' 't-value' 'p-value' 'q-value'],string(ExternalVariable_Label),'A1')
        xlswrite(filename_results,FTMS_CorrelationResults,string(ExternalVariable_Label),'A2');
    end

else
    p(p >= 1-format.CL_alpha/100) = NaN;
    SpearmanResults=[r',t',p'];
    FTMS_CorrelationResults=[FTMS_Formulas_Data, SpearmanResults];
    FTMS_CorrelationResults(any(isnan(FTMS_CorrelationResults), 2), :) = [];
    
    if ~isempty(FTMS_CorrelationResults)
        % Export data in the Excel file
        xlswrite(filename_results,[FTMS_Formulas_Labels 'r-value' 't-value' 'p-value'],string(ExternalVariable_Label),'A1')
        xlswrite(filename_results,FTMS_CorrelationResults,string(ExternalVariable_Label),'A2'); 
    end
    
end

figure
hold on
    xlim([0 1.2]); ylim([0 2.5]);
    xlabel('O/C Ratio','fontsize',13,'Position',[0.60000059604647,-0.17,-0.999999999999986]); 
    ylabel('H/C Ratio','fontsize',13,'Position',[-0.12,1.2500011920929,-0.999999999999986]);
    title(['Spearman Correlation with ' ExternalVariable_Label],'Color','k','fontsize',13,'Position',[0.600000814387664,2.55,0])

    FTMS_CorrelationResults_POS=FTMS_CorrelationResults(FTMS_CorrelationResults(:,14) >=0,:);   % Column 14 = rho value
    FTMS_CorrelationResults_NEG=FTMS_CorrelationResults(FTMS_CorrelationResults(:,14) < 0,:);   % Column 14 = rho value
    OC_ALL=FTMS_Formulas_Data(:,9);            HC_ALL=FTMS_Formulas_Data(:,10);                 % Column 9 = O/C, 10 = H/C
    OC_POS=FTMS_CorrelationResults_POS(:,9);   HC_POS=FTMS_CorrelationResults_POS(:,10);        % Column 9 = O/C, 10 = H/C
    OC_NEG=FTMS_CorrelationResults_NEG(:,9);   HC_NEG=FTMS_CorrelationResults_NEG(:,10);        % Column 9 = O/C, 10 = H/C

    String_POS = ['POS (' num2str(size(OC_POS,1)) ', ' num2str(round(size(OC_POS,1)*100/size(OC_ALL,1))) '%)' '  '];
    String_NEG = ['NEG (' num2str(size(OC_NEG,1)) ', ' num2str(round(size(OC_NEG,1)*100/size(OC_ALL,1))) '%)' '  '];
    String_ALL = ['ALL (' num2str(size(OC_ALL,1)) ', ' num2str(100) '%)' '  '];        

    A=plot(OC_NEG,HC_NEG,'b.','Markersize',11,'DisplayName',String_NEG);  % Blue Negative
    B=plot(OC_POS,HC_POS,'r.','MarkerSize',11,'DisplayName',String_POS);  % Red Positive
    C=plot(OC_ALL,HC_ALL,'y.','Markersize',13, 'DisplayName',String_ALL); % Gray all
    set(C,'Color',[0.50,0.50,0.50]);
    H=legend([A C B]);
    if ~isempty(OC_POS) && ~isempty(OC_NEG)    % Both + and - correllations exist
        h=get(gca,'Children');
        set(gca,'Children',[h(3) h(2) h(1)])
    elseif isempty(OC_POS) && ~isempty(OC_NEG) % Only negative correlation exist
        h=get(gca,'Children');
        set(gca,'Children',[h(2) h(1)])
    elseif ~isempty(OC_POS) && isempty(OC_NEG) % Only positive correlation exist
        h=get(gca,'Children');
        set(gca,'Children',[h(2) h(1)])
    else
        h=get(gca,'Children');
        set(gca,'Children',[h(1)])        
    end
    h1=refline(-1,2);               % AImod=0
    h2=refline(-0.3845,1.0584);     % AImod=0.5
    h3=refline(-0.3740,0.7551);     % AImod=0.67
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
    H.Position = [4.25003334776883,3.48916666666667,1.499999971129,0.552083318432173];
    H.Orientation = 'horizontal';
    H.NumColumns = 1;
hold off
print(gcf,['vK_Spearman_' ExternalVariable_Label],'-dpng','-r500');
end
        
disp(['Spearman correlation analysis - Complete! (' num2str(toc) ' seconds)'])        
end