function FTMS_Compare3(filename1,filename2,filename3)
%% Description 

% This code compares three samples. It finds common formulas between the 
% three samples and labels the remaining as unique formulas. It then plots:
    % 1) all common formulas on vK diagrams color-coded based on rel. magnitude
    % 2) formulas unique to 1 or common to 2 samples color-coded based on 
    %    number of C atoms (Z_scale = 'carbon') or N/C ratio (Z_scale = 'nitrogen')

% The code also exports the different groups of formulas (all, common, unique) 
% for each sample and employes FTMS_Metrics on each set. 

Z_scale='carbon'; % Pick 'carbon' or 'nitrogen'

% Example:
% FTMS_Compare3('Sample 1_Final.xlsx','Sample 2_Final.xlsx','Sample 3_Final.xlsx')

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

% Load data 
[Data1,TXT] = xlsread(filename1); titles=TXT(1,:);
Data2 = xlsread(filename2);  
Data3 = xlsread(filename3); 

Data1(:,format.Column_ExactMass)=Data1(:,format.Column_C)*format.Mass_12C+Data1(:,format.Column_H)*format.Mass_1H+Data1(:,format.Column_N)*format.Mass_14N+...
    Data1(:,format.Column_O)*format.Mass_16O+Data1(:,format.Column_S)*format.Mass_32S+Data1(:,format.Column_P)*format.Mass_31P+Data1(:,format.Column_E)*format.Mass_Heteroelement;

Data2(:,format.Column_ExactMass)=Data2(:,format.Column_C)*format.Mass_12C+Data2(:,format.Column_H)*format.Mass_1H+Data2(:,format.Column_N)*format.Mass_14N+...
    Data2(:,format.Column_O)*format.Mass_16O+Data2(:,format.Column_S)*format.Mass_32S+Data2(:,format.Column_P)*format.Mass_31P+Data2(:,format.Column_E)*format.Mass_Heteroelement;

Data3(:,format.Column_ExactMass)=Data3(:,format.Column_C)*format.Mass_12C+Data3(:,format.Column_H)*format.Mass_1H+Data3(:,format.Column_N)*format.Mass_14N+...
    Data3(:,format.Column_O)*format.Mass_16O+Data3(:,format.Column_S)*format.Mass_32S+Data3(:,format.Column_P)*format.Mass_31P+Data3(:,format.Column_E)*format.Mass_Heteroelement;

Data1(:,format.Column_ExactMass)=round(Data1(:,format.Column_ExactMass),format.Precision);
Data2(:,format.Column_ExactMass)=round(Data2(:,format.Column_ExactMass),format.Precision);
Data3(:,format.Column_ExactMass)=round(Data3(:,format.Column_ExactMass),format.Precision);

Data1 = sortrows(Data1,format.Column_ExactMass);
Data2 = sortrows(Data2,format.Column_ExactMass);
Data3 = sortrows(Data3,format.Column_ExactMass); 

RelH_Data1=Data1(:,format.Column_Magnitude)/sum(Data1(:,format.Column_Magnitude));
RelH_Data2=Data2(:,format.Column_Magnitude)/sum(Data2(:,format.Column_Magnitude));
RelH_Data3=Data3(:,format.Column_Magnitude)/sum(Data3(:,format.Column_Magnitude));

% Identify Common 
Common=intersect(Data1(:,format.Column_ExactMass),Data2(:,format.Column_ExactMass));
Common=intersect(Common,Data3(:,format.Column_ExactMass));

% Find the row location of common formulas
[~,~,index1common]=intersect(Common,Data1(:,format.Column_ExactMass)); 
[~,~,index2common]=intersect(Common,Data2(:,format.Column_ExactMass));
[~,~,index3common]=intersect(Common,Data3(:,format.Column_ExactMass));

% Check
if size(index1common,1) ~= size(index2common,1)
    error('Error! You have the same molecular formula assigned to more than one peak! Check your formula lists for ExactMass duplicates!')
end

if size(index1common,1) ~= size(index3common,1)
    error('Error! You have the same molecular formula assigned to more than one peak! Check your formula lists for ExactMass duplicates!')
end

% Extract common formulas
Common1=Data1(index1common,:);
Common2=Data2(index2common,:);
Common3=Data3(index3common,:);

% Identify Unique
index1unique=[1:1:size(Data1(:,format.Column_ExactMass),1)]'; index1unique(index1common)=[];
index2unique=[1:1:size(Data2(:,format.Column_ExactMass),1)]'; index2unique(index2common)=[];
index3unique=[1:1:size(Data3(:,format.Column_ExactMass),1)]'; index3unique(index3common)=[];

% Extract unique formulas
Unique1=Data1(index1unique,:);
Unique2=Data2(index2unique,:);
Unique3=Data3(index3unique,:);

% Extract data for vK diagrams
HC_Common1=Common1(:,format.Column_HC);  OC_Common1=Common1(:,format.Column_OC);  RelH_Common1=Common1(:,format.Column_Magnitude)/sum(Data1(:,format.Column_Magnitude));
HC_Common2=Common2(:,format.Column_HC);  OC_Common2=Common2(:,format.Column_OC);  RelH_Common2=Common2(:,format.Column_Magnitude)/sum(Data2(:,format.Column_Magnitude));
HC_Common3=Common3(:,format.Column_HC);  OC_Common3=Common3(:,format.Column_OC);  RelH_Common3=Common3(:,format.Column_Magnitude)/sum(Data3(:,format.Column_Magnitude));

HC_Unique1=Unique1(:,format.Column_HC);  OC_Unique1=Unique1(:,format.Column_OC);  RelH_Unique1=Unique1(:,format.Column_Magnitude)/sum(Data1(:,format.Column_Magnitude));
HC_Unique2=Unique2(:,format.Column_HC);  OC_Unique2=Unique2(:,format.Column_OC);  RelH_Unique2=Unique2(:,format.Column_Magnitude)/sum(Data2(:,format.Column_Magnitude));
HC_Unique3=Unique3(:,format.Column_HC);  OC_Unique3=Unique3(:,format.Column_OC);  RelH_Unique3=Unique3(:,format.Column_Magnitude)/sum(Data3(:,format.Column_Magnitude));

if strcmpi(Z_scale,'carbon')
    C_Unique1=Unique1(:,format.Column_C);
    C_Unique2=Unique2(:,format.Column_C);
    C_Unique3=Unique3(:,format.Column_C);    
elseif strcmpi(Z_scale,'nitrogen')
    NC_Unique1=Unique1(:,format.Column_N)./Unique1(:,format.Column_C);
    NC_Unique2=Unique2(:,format.Column_N)./Unique2(:,format.Column_C);
    NC_Unique3=Unique3(:,format.Column_N)./Unique3(:,format.Column_C);
else
    error('Error!Wrong argument for Z_scale!')
end

%% Figures

figurename=[char(filename1(1:end-5)) '_' char(filename2(1:end-5)) '_' char(filename3(1:end-5))];
figure(1) % Multiple vK of all common formulas
    hold on
    set(gcf,'WindowState','maximize')
    view(25,18)
    set(gca,'OuterPosition',[-0.0390625,0,1.000000000000001,1]);

    th = 0:pi/50:2*pi; 
    xunit = 0.2 * cos(th) + 0.6; 
    yunit = 0.2 * sin(th) + 1;
    zunit_Sample1 = zeros(size(xunit)); 
    zunit_Sample2 = zeros(size(xunit)) + 1;
    zunit_Sample3 = zeros(size(xunit)) + 2; 
    
    Sample1_Common=zeros(size(OC_Common1));
    Sample2_Common=(zeros(size(OC_Common1)))+1;
    Sample3_Common=(zeros(size(OC_Common1)))+2;

    plot3(zunit_Sample1,xunit,yunit,'k-','LineWidth',15);   % Circle1
    plot3(zunit_Sample2,xunit,yunit,'k-','LineWidth',15);   % Circle2
    plot3(zunit_Sample3,xunit,yunit,'k-','LineWidth',15);   % Circle3

    scatter3(Sample1_Common,OC_Common1,HC_Common1,70,log(RelH_Common1),'filled','o');
    scatter3(Sample2_Common,OC_Common2,HC_Common2,70,log(RelH_Common2),'filled','o');
    scatter3(Sample3_Common,OC_Common3,HC_Common3,70,log(RelH_Common3),'filled','o');

    color1=[0,0,0.562500000000000;0,0,0.605255670168183;0,0,0.648011340336366;0,0,0.690767066045241;0,0,0.733522740277377;0,0,0.776278410445560;0,0,0.819034080613743;0,0,0.861789750781926;0,0,0.904545469717546;0,0,0.947301150722937;0,0,0.990056820891120;0.0468749999999999,0.0468749999999999,0.953125000000000;0.107954545454545,0.107954545454545,0.892045454545455;0.169034090909091,0.169034090909091,0.830965909090909;0.230113636363636,0.230113636363636,0.769886363636364;0.291193181818182,0.291193181818182,0.708806818181819;0.352272727272727,0.352272727272727,0.647727272727273;0.413352272727273,0.413352272727273,0.586647727272728;0.474431818181818,0.474431818181818,0.525568181818182;0.535511363636363,0.535511363636363,0.464488636363637;0.596590909090909,0.596590909090909,0.403409090909091;0.657670454545454,0.657670454545454,0.342329545454546;0.718750000000000,0.718750000000000,0.281250000000000;0.779829545454545,0.779829545454545,0.220170454545455;0.840909090909090,0.840909090909090,0.159090909090910;0.901988636363636,0.901988636363636,0.0980113636363642;0.963068181818181,0.963068181818181,0.0369318181818188;1,0.972402587532998,0;1,0.902597389437936,0;1,0.832792207598687,0;1,0.762987006794324,0;1,0.693181827664376,0;1,0.623376624150711,0;1,0.553571447730065,0;1,0.483766234733843,0;1,0.413961044089362,0;1,0.344155853444881,0;1,0.274350660429761,0;1,0.204545457254758,0;1,0.134740265763620,0;1,0.0649350678378899,0;0.998517786914652,0,0;0.977272732691332,0,0;0.956521749496460,0,0;0.933695667982102,0,0;0.910869586467743,0,0;0.888043496012688,0,0;0.865217363834381,0,0;0.842391282320023,0,0;0.819565200805664,0,0;0.796739119291306,0,0;0.773913037776948,0,0;0.751086956262589,0,0;0.728260874748231,0,0;0.705434793233872,0,0;0.682608711719514,0,0;0.659782630205155,0,0;0.636956548690797,0,0;0.614130422472955,0,0;0.591304326057435,0,0;0.568478244543077,0,0;0.545652163028718,0,0;0.522826081514360,0,0;0.500000000000001,0,0];
    colormap(color1)
    
    c1=colorbar;
    caxis([min([log(RelH_Common1);log(RelH_Common2);log(RelH_Common3)]) max([log(RelH_Common1);log(RelH_Common2);log(RelH_Common3)])]);
    set(c1,'Ticks',[]);
    set(c1,'Position',[0.94,0.11,0.011111111111111,0.815]);
    set(get(c1,'title'),'string',{'Most';'Abundant'}','fontsize',20);
    set(get(c1,'XLabel'),'String',{'Least';'Abundant'},'Rotation',0,'fontsize',20,'Units','points','Position',[8,0,0],'Color','k')

    xlabel('Sample Number','Position',[1.07779506678768,-0.039091183488047,-0.089990101745784]); 
    ylabel('O/C');
    zlabel('H/C')
    title(figurename,'Interpreter','none');
    ylim([0 1.2]);      % H/C
    zlim([0.2 2.25]);   % O/C

    set(gca,'fontweight','bold','XTickLabelRotation',30)
    xticks([0 1 2]); 
    xticklabels({'Sample 1','Sample 2','Sample 3'});
    yticks([0 0.3 0.6 0.9 1.2]);
    zticks([0.5 1 1.5 2 2.5]);
    set(gca,'FontSize',15)
    set(gca,'ycolor','k')
    set(gca,'xcolor','k')
    set(gca,'LineWidth',3)
    hold off 
print(gcf,['Compare3_Common_' figurename],'-dtiffn','-r400');

%% Multiple vK of all unique formulas

figure(2) 
    hold on
    set(gcf,'WindowState','maximize')
    view(25,18)

    Line_x=0:0.2:0.8;
    Line_y_Condensed=(-0.3740*Line_x)+0.7551;

    plot3(zeros(size(Line_x)),Line_x,Line_y_Condensed,'k-','LineWidth',3);
    plot3(zeros(size(Line_x)),Line_x,zeros(size(Line_x))+0.3,'k-','LineWidth',3);
    plot3([0 0 0],[0.8 0.8 0.8],[0.3 0.4 0.4559],'k-','LineWidth',3)

    plot3(zeros(size(Line_x))+1,Line_x,Line_y_Condensed,'k-','LineWidth',3);
    plot3(zeros(size(Line_x))+1,Line_x,zeros(size(Line_x))+0.3,'k-','LineWidth',3);
    plot3([1 1 1],[0 0 0],[0.3 0.5 0.7551],'k-','LineWidth',3)
    plot3([1 1 1],[0.8 0.8 0.8],[0.3 0.4 0.4559],'k-','LineWidth',3)

    plot3(zeros(size(Line_x))+2,Line_x,Line_y_Condensed,'k-','LineWidth',3);
    plot3(zeros(size(Line_x))+2,Line_x,zeros(size(Line_x))+0.3,'k-','LineWidth',3);
    plot3([2 2 2],[0 0 0],[0.3 0.5 0.7551],'k-','LineWidth',3)
    plot3([2 2 2],[0.8 0.8 0.8],[0.3 0.4 0.4559],'k-','LineWidth',3)
    
    Sample1_Unique=zeros(size(OC_Unique1));
    Sample2_Unique=(zeros(size(OC_Unique2)))+1;
    Sample3_Unique=(zeros(size(OC_Unique3)))+2;

    if strcmpi(Z_scale,'nitrogen')
        scatter3(Sample1_Unique,OC_Unique1,HC_Unique1,30,NC_Unique1,'filled','o');
        scatter3(Sample2_Unique,OC_Unique2,HC_Unique2,30,NC_Unique2,'filled','o');
        scatter3(Sample3_Unique,OC_Unique3,HC_Unique3,30,NC_Unique3,'filled','o');
        colorscheme=[0,0,0;0,0,1;0,0,1;0.392156869173050,0.831372559070587,0.0745098069310188;0.392156869173050,0.831372559070587,0.0745098069310188;0.392156869173050,0.831372559070587,0.0745098069310188;1,0,0;1,0,0;1,0,0;1,0,0];
        colormap(colorscheme);
        c2=colorbar;
        caxis([min([NC_Unique1;NC_Unique2;NC_Unique3]) max([NC_Unique1;NC_Unique2;NC_Unique3])]);
        set(c2,'Position',[0.95,0.148481714893304,0.011111111111111,0.776518285106695]);
        set(get(c2,'title'),'string', 'N/C Ratio','fontsize',20,'Position',[11,590,0]);
    elseif strcmpi(Z_scale,'carbon')
        scatter3(Sample1_Unique,OC_Unique1,HC_Unique1,30,C_Unique1,'filled','o');
        scatter3(Sample2_Unique,OC_Unique2,HC_Unique2,30,C_Unique2,'filled','o');
        scatter3(Sample3_Unique,OC_Unique3,HC_Unique3,30,C_Unique3,'filled','o');
        colorscheme=[0,0,0;0,0,1;0,0,1;0.392156869173050,0.831372559070587,0.0745098069310188;0.392156869173050,0.831372559070587,0.0745098069310188;0.392156869173050,0.831372559070587,0.0745098069310188;1,0,0;1,0,0;1,0,0;1,0,0];
        colormap(colorscheme);
        c2=colorbar;
        caxis([0 60]);
        set(c2,'Ticks',[6,12,18,24,30,36,42,48,54,60]);
        set(c2,'Position',[0.96,0.148481714893304,0.011111111111111,0.776518285106695]);
        set(get(c2,'title'),'string', 'C number','fontsize',20,'Position',[8,580,0]);
    else
        error('Error!Wrong argument for Z_scale!')
    end

    set(c2,'FontSize',20);
    xlabel('Sample Number','Position',[1.07779506678768,-0.039091183488047,-0.089990101745784]); 
    ylabel('O/C');
    zlabel('H/C');
    title(figurename,'Interpreter','none');
    ylim([0 1.2]);     % H/C
    zlim([0.2 2.25]);  % O/C
    set(gca,'fontweight','bold','XTickLabelRotation',30)
    xticks([0 1 2]); 
    xticklabels({'Sample 1','Sample 2','Sample 3'});
    yticks([0 0.3 0.6 0.9 1.2]);
    zticks([0.5 1 1.5 2 2.5]);
    set(gca,'FontSize',15)
    set(gca,'ycolor','k')
    set(gca,'xcolor','k')
    set(gca,'LineWidth',3)
    hold off
if strcmpi(Z_scale,'carbon')
    print(gcf,['Compare3_Unique_C_' figurename],'-dtiffn','-r400');
elseif strcmpi(Z_scale,'nitrogen')
    print(gcf,['Compare3_Unique_N_' figurename],'-dtiffn','-r400');
else
    error('Error!Wrong argument for Z_scale!')
end
    

%% Data Export 

warning('off','MATLAB:xlswrite:AddSheet');
exportname=['Comparison3_' char(figurename) '.xlsx'];

titles=[titles 'RelMagn'];

xlswrite(exportname,titles,'Data1','A1');   xlswrite(exportname,[Data1,RelH_Data1],'Data1','A2');
xlswrite(exportname,titles,'Data2','A1');   xlswrite(exportname,[Data2,RelH_Data2],'Data2','A2');
xlswrite(exportname,titles,'Data3','A1');   xlswrite(exportname,[Data3,RelH_Data3],'Data3','A2');

xlswrite(exportname,titles,'Common1','A1'); xlswrite(exportname,[Common1,RelH_Common1],'Common1','A2');
xlswrite(exportname,titles,'Common2','A1'); xlswrite(exportname,[Common2,RelH_Common2],'Common2','A2');
xlswrite(exportname,titles,'Common3','A1'); xlswrite(exportname,[Common3,RelH_Common3],'Common3','A2');

xlswrite(exportname,titles,'Unique1','A1'); xlswrite(exportname,[Unique1,RelH_Unique1],'Unique1','A2');
xlswrite(exportname,titles,'Unique2','A1'); xlswrite(exportname,[Unique2,RelH_Unique2],'Unique2','A2');
xlswrite(exportname,titles,'Unique3','A1'); xlswrite(exportname,[Unique3,RelH_Unique3],'Unique3','A2');

loc1='W3'; loc2='W7';
FTMS_Metrics(exportname,'Data1',loc1,loc2,true)
FTMS_Metrics(exportname,'Data2',loc1,loc2,true)
FTMS_Metrics(exportname,'Data3',loc1,loc2,true)

FTMS_Metrics(exportname,'Common1',loc1,loc2,true)
FTMS_Metrics(exportname,'Common2',loc1,loc2,true)
FTMS_Metrics(exportname,'Common3',loc1,loc2,true)

FTMS_Metrics(exportname,'Unique1',loc1,loc2,true)
FTMS_Metrics(exportname,'Unique2',loc1,loc2,true)
FTMS_Metrics(exportname,'Unique3',loc1,loc2,true)

disp(['Finished comparing ' char(filename1) ', ' char(filename2) ' and ' char(filename3)])
end