function FTMS_Compare(filename1,filename2)
%% Description

% This code compares two samples. It finds common formulas between the 
% two samples and finds the unique formulas of each. It then plots:
    % 1) all formulas and unique formulas on vK and H/C vs ExactMass plots
    % 2) 3D Magnitude vK plots of common formulas of both samples as well as percent difference
    % 3) Venn diagram with statistics as well as a Magnitude Comparison per Sleighter et al., 2012. 

% The code also exports the different groups of formulas (all, common, unique) 
% for each sample and employes FTMS_Metrics on each set. 

% Example:
% FTMS_Compare('Sample 1_Final.xlsx','Sample 2_Final.xlsx')

%% Copyright and License Notices: 

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
[Data1,TXT] = xlsread(filename1); titles=TXT(1,:);
Data2 = xlsread(filename2); 

Data1(:,format.Column_ExactMass)=Data1(:,format.Column_C)*format.Mass_12C+Data1(:,format.Column_H)*format.Mass_1H+Data1(:,format.Column_N)*format.Mass_14N+...
    Data1(:,format.Column_O)*format.Mass_16O+Data1(:,format.Column_S)*format.Mass_32S+Data1(:,format.Column_P)*format.Mass_31P+Data1(:,format.Column_E)*format.Mass_Heteroelement;

Data2(:,format.Column_ExactMass)=Data2(:,format.Column_C)*format.Mass_12C+Data2(:,format.Column_H)*format.Mass_1H+Data2(:,format.Column_N)*format.Mass_14N+...
    Data2(:,format.Column_O)*format.Mass_16O+Data2(:,format.Column_S)*format.Mass_32S+Data2(:,format.Column_P)*format.Mass_31P+Data2(:,format.Column_E)*format.Mass_Heteroelement;

Data1(:,format.Column_ExactMass)=round(Data1(:,format.Column_ExactMass),format.Precision);
Data2(:,format.Column_ExactMass)=round(Data2(:,format.Column_ExactMass),format.Precision);

Data1 = sortrows(Data1,format.Column_ExactMass); 
Data2 = sortrows(Data2,format.Column_ExactMass); 

RelH_Data1=Data1(:,format.Column_Magnitude)/sum(Data1(:,format.Column_Magnitude));
RelH_Data2=Data2(:,format.Column_Magnitude)/sum(Data2(:,format.Column_Magnitude));

%Identify Common 
Common=intersect(Data1(:,format.Column_ExactMass),Data2(:,format.Column_ExactMass));

%Find the row location of common formulas
[~,~,index1common]=intersect(Common,Data1(:,format.Column_ExactMass)); 
[~,~,index2common]=intersect(Common,Data2(:,format.Column_ExactMass));

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

% Extract data for vK diagrams
HC_Data1=Data1(:,format.Column_HC);      OC_Data1=Data1(:,format.Column_OC);
HC_Data2=Data2(:,format.Column_HC);      OC_Data2=Data2(:,format.Column_OC);

HC_Common1=Common1(:,format.Column_HC);  OC_Common1=Common1(:,format.Column_OC);  RelH_Common1=Common1(:,format.Column_Magnitude)/sum(Data1(:,format.Column_Magnitude));
HC_Common2=Common2(:,format.Column_HC);  OC_Common2=Common2(:,format.Column_OC);  RelH_Common2=Common2(:,format.Column_Magnitude)/sum(Data2(:,format.Column_Magnitude));

HC_Unique1=Unique1(:,format.Column_HC);  OC_Unique1=Unique1(:,format.Column_OC);  RelH_Unique1=Unique1(:,format.Column_Magnitude)/sum(Data1(:,format.Column_Magnitude));
HC_Unique2=Unique2(:,format.Column_HC);  OC_Unique2=Unique2(:,format.Column_OC);  RelH_Unique2=Unique2(:,format.Column_Magnitude)/sum(Data2(:,format.Column_Magnitude));

%% Statistics

N_Data1=size(Data1,1);      N_Data2=size(Data2,1);     
N_Common1=size(Common1,1);  N_Common2=size(Common2,1);  
N_Unique1=size(Unique1,1);  N_Unique2=size(Unique2,1);  

%% Figures

figurename=[char(filename1(1:end-5)) '_' char(filename2(1:end-5))];

figure(1)
set(gcf,'WindowState','maximized')
    subplot(2,2,1) % vK stacked two samples
        hold on
        title('vK All Formulas','Color','k','fontsize',13,'Position',[0.600000711628944,2.5669,1.4e-14]);
        xlabel('O/C','fontsize',13); 
        ylabel('H/C','fontsize',13);
        xlim([0 1.2]); ylim([0 2.5]);
        
        String_A=[char(filename1(1:end-5)) ' (' num2str(N_Data1) ' total formulas)'];
        String_B=[char(filename2(1:end-5)) ' (' num2str(N_Data2) ' total formulas)'];
        
        A = plot(OC_Data1,HC_Data1,'ro','MarkerSize',6,'DisplayName',String_A);
        B = plot(OC_Data2,HC_Data2,'b.','MarkerSize',7,'DisplayName',String_B);
        
        l=legend([A B]);
        set(l,'Interpreter', 'none');
        set(l,'FontWeight', 'bold');
        set(l,'FontSize', 12);
        
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
        set(gca,'FontSize',15)
        set(gca,'OuterPosition',[0.07896013864818,0.51945507182499,0.392990631997137,0.475861896967504]);
        hold off
    subplot(2,2,2) % vK unique only
        hold on
        title('vK Unique Formulas','Color','k','fontweight','bold','fontsize',13,'Position',[0.600002160668396,2.5669,1.4e-14]);
        xlabel('O/C','fontsize',13); 
        ylabel('H/C','fontsize',13);
        xlim([0 1.2]); ylim([0 2.5]);
        
        String_C=[char(filename1(1:end-5)) ' (' num2str(N_Unique1) ' unique, ' num2str(N_Common1) ' common formulas)'];
        String_D=[char(filename2(1:end-5)) ' (' num2str(N_Unique2) ' unique, ' num2str(N_Common2) ' common formulas)'];
        
        C = plot(OC_Unique1,HC_Unique1,'ro','MarkerSize',6,'DisplayName',String_C);
        D = plot(OC_Unique2,HC_Unique2,'b.','MarkerSize',7,'DisplayName',String_D);
        
        l=legend([C D]);
        set(l,'Interpreter', 'none');
        set(l,'FontWeight', 'bold');
        set(l,'FontSize', 12);
        
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
        set(gca,'FontSize',15)
        set(gca,'OuterPosition',[0.534263458060926,0.51945507182499,0.409653637501739,0.475861896967504]);
        hold off
    subplot(2,2,3) % H/C vs MW stacked two samples
        hold on
        title('H/C vs MW All Formulas','Color','k','fontsize',13);
        xlabel('Molecular Weight','fontsize',13); 
        ylabel('H/C','fontsize',13);
        xlim([150 1000]); ylim([0 2.5]);
        xticks([200 400 600 800 1000]);
        
        plot(Data1(:,format.Column_ExactMass),HC_Data1,'ro','MarkerSize',6);
        plot(Data2(:,format.Column_ExactMass),HC_Data2,'b.','MarkerSize',7);
        set(gca,'fontweight','bold')
        set(gca,'ycolor','k')
        set(gca,'xcolor','k')
        set(gca,'OuterPosition',[0.07896013864818,0.009345796063682,0.397417701911665,0.476604117339637]);
        set(gca,'FontSize',15)
        hold off
    subplot(2,2,4) % H/C vs MW stacked two samples
        hold on
        title('H/C vs MW Unique Formulas','Color','k','fontsize',13);
        xlabel('Molecular Weight','fontsize',13); 
        ylabel('H/C','fontsize',13);
        xlim([150 1000]); ylim([0 2.5]);
        xticks([200 400 600 800 1000]);
        
        plot(Unique1(:,format.Column_ExactMass),HC_Unique1,'ro','MarkerSize',6);
        plot(Unique2(:,format.Column_ExactMass),HC_Unique2,'b.','MarkerSize',7);
        set(gca,'fontweight','bold')
        set(gca,'ycolor','k')
        set(gca,'xcolor','k')
        set(gca,'FontSize',15)
        set(gca,'OuterPosition',[0.534263458060926,0.009345796063682,0.409653637501739,0.476604117339637]);
        hold off        
print(gcf,['Compare_Unique_' figurename],'-dpng','-r300');

figure(2) % vK - 3D plots of common formulas
set(gcf,'WindowState','maximize') 
    ax(1)=subplot(1,3,1); % Vk 3D Sample 1
        hold on
        title(char(filename1(1:end-5)),'Interpreter','none','fontsize',13,'Position',[0.600001048914274,2.55,-5.1633191108703]); 

        IntLog_Common1=log(RelH_Common1);
        IntLog_Common2=log(RelH_Common2);
        bottom=min(min(min(IntLog_Common1)),min(min(IntLog_Common2)));
        top=max(max(max(IntLog_Common1)),max(max(IntLog_Common2)));

        scatter3(OC_Common1(:),HC_Common1(:),IntLog_Common1(:),30,IntLog_Common1(:),'filled','o');
        view(0,90) %(45,30)

        c1=colorbar;
        set(get(c1,'title'),'string', 'log(RelMagn)','fontsize',10);
        caxis manual
        caxis([bottom top]);

        r=circle3(0.6,1,0,0.2);
        r.LineWidth=3;

        grid on        
        originalSize1 = get(gca, 'Position');
        xlabel('O/C Ratio','fontsize',13,'Position',[0.600000596046471,-0.16,-10.3266382217407]); 
        ylabel('H/C Ratio','fontsize',13,'Position',[-0.2,1.2500011920929,-10.3266382217407]); 
        xlim([0 1.2]); ylim([0 2.5]);
        xticks([0.3 0.6 0.9 1.2])
        set(gca,'fontweight','bold')
        set(gca,'ycolor','k')
        set(gca,'xcolor','k')
        axis manual
        set(gca, 'Position', originalSize1); 
        hold off
    ax(2)=subplot(1,3,2); % Vk 3D Sample 2
        hold on
        title(char(filename2(1:end-5)),'Interpreter','none','fontsize',13,'Position',[0.600001048914274,2.55,-5.1633191108703]);

        scatter3(OC_Common1(:),HC_Common1(:),IntLog_Common2(:),30,IntLog_Common2(:),'filled','o');
        view(0,90)

        c2=colorbar;
        set(get(c2,'title'),'string', 'log(RelMagn)','fontsize',10);
        caxis manual
        caxis([bottom top]);

        r=circle3(0.6,1,0,0.2);
        r.LineWidth=3;

        grid on   
        originalSize2 = get(gca, 'Position');
        xlabel('O/C Ratio','fontsize',13,'Position',[0.6,-0.16,-5.1866])
        %ylabel('H/C Ratio','fontsize',13) % No y-label for cleaner figure
        xlim([0 1.2]); ylim([0 2.5]);
        xticks([0.3 0.6 0.9 1.2]);
        set(gca,'fontweight','bold');
        set(gca,'ycolor','k');
        set(gca,'xcolor','k');
        axis manual
        set(gca, 'Position', originalSize2);
        hold off    
    ax(3)=subplot(1,3,3); % Vk 3D Difference
        hold on
        title('% Change','Interpreter','none','fontsize',13,'Position',[0.600001048914274,2.55,-5.1633191108703]);

        RelInt_Diff=((RelH_Common2-RelH_Common1)./RelH_Common1)*100;

        scatter3(OC_Common1(:),HC_Common1(:),RelInt_Diff(:),30,RelInt_Diff(:),'filled','o');
        view(0,90)

        c3=colorbar;
        set(get(c3,'title'),'string', '%Diff','fontsize',10);
        caxis manual
        caxis([-100 400]);

        r=circle3(0.6,1,600,0.2);
        r.LineWidth=3;        

        grid on  
        originalSize3 = get(gca, 'Position');
        xlabel('O/C Ratio','fontsize',13,'Position',[0.6,-0.16,-5.1866]) 
        %ylabel('H/C Ratio','fontsize',13) % No y-label for cleaner figure
        xlim([0 1.2]); ylim([0 2.5]);
        xticks([0.3 0.6 0.9 1.2]);
        set(gca,'fontweight','bold')
        set(gca,'ycolor','k')
        set(gca,'xcolor','k')
        axis manual
        set(gca, 'Position', originalSize3); 
        hold off    
c=[0,0,1;0.00573152350261807,0.0310708899050951,0.993061840534210;0.0114630470052361,0.0621417798101902,0.986123681068420;0.0171945709735155,0.0932126715779305,0.979185521602631;0.0229260940104723,0.124283559620380,0.972247362136841;0.0286576189100742,0.155354455113411,0.965309202671051;0.0343891419470310,0.186425343155861,0.958371043205261;0.0401206649839878,0.217496231198311,0.951432883739471;0.0458521880209446,0.248567119240761,0.944494724273682;0.0515837110579014,0.279638022184372,0.937556564807892;0.0573152378201485,0.310708910226822,0.930618405342102;0.0630467608571053,0.341779798269272,0.923680245876312;0.0687782838940620,0.372850686311722,0.916742086410523;0.0745098069310188,0.403921574354172,0.909803926944733;0.0802413299679756,0.434992462396622,0.902865767478943;0.0859728530049324,0.466063350439072,0.895927608013153;0.0917043760418892,0.497134238481522,0.888989448547363;0.0974358990788460,0.528205156326294,0.882051289081574;0.103167422115803,0.559276044368744,0.875113129615784;0.108898945152760,0.590346932411194,0.868174970149994;0.114630475640297,0.621417820453644,0.861236810684204;0.120361998677254,0.652488708496094,0.854298651218414;0.126093521714211,0.683559596538544,0.847360491752625;0.131825044751167,0.714630484580994,0.840422332286835;0.137556567788124,0.745701372623444,0.833484172821045;0.143288090825081,0.776772260665894,0.826546013355255;0.149019613862038,0.807843148708344,0.819607853889465;0.138375356793404,0.821568608283997,0.832493007183075;0.127731099724770,0.835294127464294,0.845378160476685;0.117086842656136,0.849019646644592,0.858263313770294;0.106442578136921,0.862745106220245,0.871148467063904;0.0957983210682869,0.876470565795898,0.884033620357513;0.0851540639996529,0.890196084976196,0.896918773651123;0.0745098069310188,0.903921604156494,0.909803926944733;0.0638655498623848,0.917647063732147,0.922689080238342;0.0532212890684605,0.931372523307800,0.935574233531952;0.0425770319998264,0.945098042488098,0.948459386825562;0.0319327749311924,0.958823561668396,0.961344540119171;0.0212885159999132,0.972549021244049,0.974229693412781;0.0106442579999566,0.986274480819702,0.987114846706390;0,1,1;0,1,0.909090936183929;0,1,0.818181812763214;0,1,0.727272748947144;0,1,0.636363625526428;0,1,0.545454561710358;0,1,0.454545468091965;0,1,0.363636374473572;0,1,0.272727280855179;0,1,0.181818187236786;0,1,0.0909090936183929;0,1,0;0.0769230797886848,1,0;0.153846159577370,1,0;0.230769231915474,1,0;0.307692319154739,1,0;0.384615391492844,1,0;0.461538463830948,1,0;0.538461565971375,1,0;0.615384638309479,1,0;0.692307710647583,1,0;0.769230782985687,1,0;0.846153855323792,1,0;0.923076927661896,1,0;1,1,0;1,0.999411761760712,0;1,0.998823523521423,0;1,0.998235285282135,0;1,0.997647047042847,0;1,0.997058808803558,0;1,0.996470570564270,0;1,0.995882332324982,0;1,0.995294094085693,0;1,0.994705855846405,0;1,0.994117617607117,0;1,0.993529379367828,0;1,0.992941141128540,0;1,0.992352962493897,0;1,0.991764724254608,0;1,0.991176486015320,0;1,0.990588247776032,0;1,0.990000009536743,0;1,0.989411771297455,0;1,0.988823533058167,0;1,0.988235294818878,0;1,0.987647056579590,0;1,0.987058818340302,0;1,0.986470580101013,0;1,0.985882341861725,0;1,0.985294103622437,0;1,0.976884782314301,0;1,0.968475401401520,0;1,0.960066080093384,0;1,0.951656758785248,0;1,0.943247437477112,0;1,0.934838056564331,0;1,0.926428735256195,0;1,0.918019413948059,0;1,0.909610033035278,0;1,0.901200711727142,0;1,0.892791390419006,0;1,0.884382069110870,0;1,0.875972688198090,0;1,0.867563366889954,0;1,0.859154045581818,0;1,0.850744664669037,0;1,0.842335343360901,0;1,0.833926022052765,0;1,0.825516700744629,0;1,0.817107319831848,0;1,0.808697998523712,0;1,0.800288677215576,0;1,0.791879296302795,0;1,0.783469974994659,0;1,0.775060653686523,0;1,0.766651272773743,0;1,0.758241951465607,0;1,0.749832630157471,0;1,0.741423308849335,0;1,0.733013927936554,0;1,0.724604606628418,0;1,0.716195285320282,0;1,0.707785904407501,0;1,0.699376583099365,0;1,0.690967261791229,0;1,0.682557940483093,0;1,0.674148559570313,0;1,0.665739238262177,0;1,0.653178095817566,0;1,0.640617012977600,0;1,0.628055870532990,0;1,0.615494787693024,0;1,0.602933645248413,0;1,0.590372502803803,0;1,0.577811419963837,0;1,0.565250277519226,0;1,0.552689194679260,0;1,0.540128052234650,0;1,0.527566969394684,0;1,0.515005826950073,0;1,0.502444684505463,0;1,0.489883601665497,0;1,0.477322459220886,0;1,0.464761346578598,0;1,0.452200233936310,0;1,0.439639121294022,0;1,0.427078008651733,0;1,0.414516896009445,0;1,0.401955753564835,0;1,0.389394640922546,0;1,0.376833528280258,0;1,0.364272415637970,0;1,0.351711302995682,0;1,0.339150190353394,0;1,0.326589047908783,0;1,0.314027935266495,0;1,0.301466822624207,0;1,0.288905709981918,0;1,0.276344597339630,0;1,0.263783484697342,0;1,0.251222342252731,0;1,0.238661229610443,0;1,0.226100116968155,0;1,0.213539004325867,0;1,0.200977876782417,0;1,0.188416764140129,0;1,0.175855651497841,0;1,0.163294523954391,0;1,0.150733411312103,0;1,0.138172298669815,0;1,0.125611171126366,0;1,0.113050058484077,0;1,0.100488938391209,0;1,0.0879278257489204,0;1,0.0753667056560516,0;1,0.0628055855631828,0;1,0.0502444691956043,0;1,0.0376833528280258,0;1,0.0251222345978022,0;1,0.0125611172989011,0;1,0,0;0.993333339691162,0,0;0.986666679382324,0,0;0.980000019073486,0,0;0.973333358764648,0,0;0.966666638851166,0,0;0.959999978542328,0,0;0.953333318233490,0,0;0.946666657924652,0,0;0.939999997615814,0,0;0.933333337306976,0,0;0.926666676998138,0,0;0.920000016689301,0,0;0.913333356380463,0,0;0.906666696071625,0,0;0.899999976158142,0,0;0.893333315849304,0,0;0.886666655540466,0,0;0.879999995231628,0,0;0.873333334922791,0,0;0.866666674613953,0,0;0.860000014305115,0,0;0.853333353996277,0,0;0.846666693687439,0,0;0.839999973773956,0,0;0.833333313465118,0,0;0.826666653156281,0,0;0.819999992847443,0,0;0.813333332538605,0,0;0.806666672229767,0,0;0.800000011920929,0,0;0.793333351612091,0,0;0.786666691303253,0,0;0.779999971389771,0,0;0.773333311080933,0,0;0.766666650772095,0,0;0.759999990463257,0,0;0.753333330154419,0,0;0.746666669845581,0,0;0.740000009536743,0,0;0.733333349227905,0,0;0.726666688919067,0,0;0.720000028610230,0,0;0.713333308696747,0,0;0.706666648387909,0,0;0.699999988079071,0,0;0.693333327770233,0,0;0.686666667461395,0,0;0.680000007152557,0,0;0.673333346843720,0,0;0.666666686534882,0,0;0.660000026226044,0,0;0.653333306312561,0,0;0.646666646003723,0,0;0.639999985694885,0,0;0.633333325386047,0,0;0.626666665077210,0,0;0.620000004768372,0,0;0.613333344459534,0,0;0.606666684150696,0,0;0.600000023841858,0,0;0.593333303928375,0,0;0.586666643619537,0,0;0.579999983310700,0,0;0.573333323001862,0,0;0.566666662693024,0,0;0.560000002384186,0,0;0.553333342075348,0,0;0.546666681766510,0,0;0.540000021457672,0,0;0.533333361148834,0,0;0.526666641235352,0,0;0.519999980926514,0,0;0.513333320617676,0,0;0.506666660308838,0,0;0.500000000000000,0,0];
colormap(ax(1),jet)
colormap(ax(2),jet)
colormap(ax(3),c)
print(gcf,['Compare_Common_' figurename],'-dpng','-r300');   

figure(3) 
set(gcf,'WindowState','maximized')
    subplot(1,2,1) % Venn diagram
        hold on
        title('Venn Diagram');
        set(gca,'OuterPosition',[-0.03605105430387,0,0.50968759720926,1]);
        xlim([-17 17]); ylim([-10 10]);
        circle(-5,0,10);
        circle(5,0,10);

        mid={['TOTAL']
            [num2str(N_Data1+N_Data2-N_Common1) ' formulas']
            [' ']
            ['COMMON']
            [num2str(N_Common1) ' formulas']
            [num2str(round(N_Common1*100/(N_Data1+N_Data2-N_Common1))) ' num%']
            [' ']
            [' ']
            [char(filename1(1:end-5))]
            [num2str(round(sum(RelH_Common1)*100)) ' magn%']
            [' ']        
            [char(filename2(1:end-5))]
            [num2str(round(sum(RelH_Common2)*100)) ' magn%']};
        
        left={[char(filename1(1:end-5))]
            [num2str(N_Data1) ' formulas']
            [' ']
            ['UNIQUE']
            [num2str(N_Unique1) ' formulas']
            [num2str(round(N_Unique1*100/(N_Data1+N_Data2-N_Common1))) ' num%']
            [' ']
            [' ']
            [' ']
            [num2str(round(sum(RelH_Unique1)*100)) ' magn%']
            [' ']
            [' ']
            [' ']};

        right={[char(filename2(1:end-5))]
            [num2str(N_Data2) ' formulas']
            [' ']
            ['UNIQUE']
            [num2str(N_Unique2) ' formulas']
            [num2str(round(N_Unique2*100/(N_Data1+N_Data2-N_Common1))) ' num%'] 
            [' ']
            [' ']
            [' ']
            [' ']
            [' ']
            [' ']
            [num2str(round(sum(RelH_Unique2)*100)) ' magn%']}; 
        
        text(0,0,mid,'Interpreter', 'none','FontSize',14,'HorizontalAlignment','center');
        text(-10,0,left,'Interpreter', 'none','FontSize',14,'HorizontalAlignment','center');
        text(10,0,right,'Interpreter', 'none','FontSize',14,'HorizontalAlignment','center');
        axis off
        hold off
    subplot(1,2,2) % Magnitude comparison
        hold on
        xlabel(char(filename1(1:end-5)),'Interpreter', 'none','fontsize',13);
        ylabel(char(filename2(1:end-5)),'Interpreter', 'none','fontsize',13);
        plot(RelH_Common1,RelH_Common2,'Marker','o','Color',[0.0745098039215686 0.623529411764706 1],'LineStyle','none');
        linefit=polyfit(RelH_Common1,RelH_Common2,1);
        R=corrcoef(RelH_Common1,RelH_Common2);
        Rsq = R(1,2).^2;
        title(['Magnitude Comparison of Common Formulas (R^2 = ' num2str(round(Rsq,3)) ')']);
        plot([min([xlim ylim]) max([xlim ylim])], [min([xlim ylim]) max([xlim ylim])], '-r','LineWidth',3); % 1:1 line
        h=refline(linefit(1),linefit(2));
        set(h,'LineWidth',3,'Color','k','LineStyle','--') % Linear Fit
        legend({'Common Formulas','One-to-One Line','Linear Fit'},'Location','northwest')
        set(gca,'fontweight','bold')
        set(gca,'ycolor','k')
        set(gca,'xcolor','k')
        hold off
print(gcf,['Compare_Stats_' figurename],'-dpng','-r300');

%% Data Export

warning('off','MATLAB:xlswrite:AddSheet');
exportname=['Comparison_' char(figurename) '.xlsx'];

titles=[titles 'RelMagn'];

xlswrite(exportname,titles,'Data1','A1');   xlswrite(exportname,[Data1,RelH_Data1],'Data1','A2');
xlswrite(exportname,titles,'Data2','A1');   xlswrite(exportname,[Data2,RelH_Data2],'Data2','A2');    

xlswrite(exportname,titles,'Common1','A1'); xlswrite(exportname,[Common1,RelH_Common1],'Common1','A2');
xlswrite(exportname,titles,'Common2','A1'); xlswrite(exportname,[Common2,RelH_Common2],'Common2','A2');

xlswrite(exportname,titles,'Unique1','A1'); xlswrite(exportname,[Unique1,RelH_Unique1],'Unique1','A2');
xlswrite(exportname,titles,'Unique2','A1'); xlswrite(exportname,[Unique2,RelH_Unique2],'Unique2','A2');

loc1='W3'; loc2='W7';
FTMS_Metrics(exportname,'Data1',loc1,loc2,true)
FTMS_Metrics(exportname,'Data2',loc1,loc2,true)

FTMS_Metrics(exportname,'Common1',loc1,loc2,true)
FTMS_Metrics(exportname,'Common2',loc1,loc2,true)

FTMS_Metrics(exportname,'Unique1',loc1,loc2,true)
FTMS_Metrics(exportname,'Unique2',loc1,loc2,true)

Table_Rows={'Total Formulas', 'Common #', 'Unique', 'Common num%', 'Unique num%','Common magn%', 'Unique magn%'};
Table_Values1=[N_Data1,N_Common1,N_Unique1,round(N_Common1*100/(N_Data1+N_Data2-N_Common1),2),round(N_Unique1*100/(N_Data1+N_Data2-N_Common1),2),round(sum(RelH_Common1)*100,2),round(sum(RelH_Unique1)*100,2)];
Table_Values2=[N_Data2,N_Common2,N_Unique2,round(N_Common2*100/(N_Data1+N_Data2-N_Common1),2),round(N_Unique2*100/(N_Data1+N_Data2-N_Common1),2),round(sum(RelH_Common2)*100,2),round(sum(RelH_Unique2)*100,2)];
Table_Values1_str=string(Table_Values1); Table_Values2_str=string(Table_Values2);

output_a=[' ',Table_Rows;'Sample 1',Table_Values1_str;'Sample 2',Table_Values2_str;];
output_b=['Slope:',string(linefit(1)),'Intercept',string(linefit(2)),' ', ' ',' ', ' '];
output_c=['R2',string(Rsq),' ', ' ',' ', ' ',' ', ' ',];
output=[output_a;output_b;output_c];

xlswrite(exportname,output);

disp(['Finished comparing ' char(filename1) ' and ' char(filename2)])
end

%% Internal functions

% Function for plotting an eclipse 
function h = eclipse(x1,x2,y1,y2,eccentricity)
    hold on
    a = (1/2) * sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2);
    b = a * sqrt(1-eccentricity^2);
    t = linspace(0, 2 * pi, 300);
    X = a * cos(t);
    Y = b * sin(t);
    angles = atan2(y2 - y1, x2 - x1);
    xunit = (x1 + x2) / 2 + X * cos(angles) - Y * sin(angles);
    yunit = (y1 + y2) / 2 + X * sin(angles) + Y * cos(angles);
    h = plot(xunit, yunit,'k-','LineWidth',2);
    hold off  
end

% Function for plotting a circle on a 2D plot
function h = circle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'k-','LineWidth',2);
    hold off    
end    

% Function for plotting a circle on a 3D plot
function h = circle3(x,y,z,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    zunit=transpose(zeros(size(th,2),1) + z);
    h = plot3(xunit, yunit,zunit,'k-','LineWidth',2);
    hold off    
end  