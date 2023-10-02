function varargout = FTMS_CompoundClass(filename)
%% Description 

% This code takes a formula list, determines the compound class for each formula
% (lignin, tannin, lipid, etc...) and type of formula (CHO, CHON, etc...)
% and does a variety of statistics (eg. % CHO, % Carbs, etc...). Also
% separates all different compound classes and formula types in different
% sheets in the produced Excel file. 

% The FTMS_Metrics code is incorporated in here and will calculate various
% metrics for each of the sheets in the "_CompoundClass.xlsx" file. 

% Example:
% FTMS_CompoundClass('Sample1_Final.xlsx')

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

[Data,TXT]=xlsread(filename);

Data=sortrows(Data,format.Column_ExactMass);
titles=TXT(1,:);

mz=Data(:,format.Column_mz);
Int=Data(:,format.Column_Magnitude);
SN=Data(:,format.Column_SN);
EstimC=Data(:,format.Column_EstimC);
C=Data(:,format.Column_C);
H=Data(:,format.Column_H);
O=Data(:,format.Column_O);
N=Data(:,format.Column_N);
S=Data(:,format.Column_S);
P=Data(:,format.Column_P);
E=Data(:,format.Column_E);
Mass=Data(:,format.Column_ExactMass);
e=Data(:,format.Column_error);
OC=Data(:,format.Column_OC);
HC=Data(:,format.Column_HC);
DBE=Data(:,format.Column_DBE);
DBEC=Data(:,format.Column_DBEC);
AImod=Data(:,format.Column_AImod);

RelH=Int./sum(Int);

% Compound Class caterogization
[r,~]=size(Mass);
w = cell(r,1);
    for x=1:r
        if E(x)>0
            w{x,1}='E'; %E-containing compounds
        elseif  AImod(x) >= 0.67 && C(x) >= 15
            w{x,1} = 'ConAC'; %Condensed aromatic compounds      
        elseif AImod(x) >= 0.67 && C(x) < 15 
            w{x,1} = 'SCA'; % Small condensed aromatics                
        elseif AImod(x) < 0.67 && HC(x) < 1.5 && OC(x) < 0.67 && OC(x) > 0.1
            w{x,1} = 'Lignin';
        elseif  HC(x)< 1.5 && OC(x)>= 0.67 && AImod(x) <0.67
            w{x,1} = 'Tannin';
        elseif HC(x)>=2 && OC(x)<0.67 && S(x)>0 
            w{x,1} = 'SA'; % Sulfonic Acids
        elseif HC(x) >= 1.5 && OC(x) >=0.55 && N(x)>0
            w{x,1} = 'AS'; % Amino Sugars            
        elseif HC(x) >= 1.5 && OC(x) >=0.67
            w{x,1} = 'Sugar';                 
        elseif HC(x)>=1.5  && OC(x)<0.55 && N(x)>0 
            w{x,1} = 'Protein';                 
        elseif HC(x)>=2 && OC(x)<0.67 && N(x)==0
            w{x,1} = 'FA'; %Fatty acids            
        elseif HC(x)<2 && HC(x)>=1.5 && OC(x) <0.67 && N(x)==0
            w{x,1} = 'Lipid';  % Includes saturated aliphatics                
        elseif HC(x)<1.5 && AImod(x)>0 && AImod(x)<0.67 && OC(x) <0.1
            w{x,1} = 'Unsaturates'; % Unsaturated compounds
        else
            w{x,1} = 'Extra';
        end
    end
    
%Formula type caterogization
[r,~]=size(Mass);    
t = cell(r,1);
    for x=1:r
        if E(x)>0
            t{x,1}='E';        
        elseif  C(x)>0 && H(x)>0 && O(x)==0 && N(x)==0 && S(x)==0 && P(x)==0 && E(x)==0
            t{x,1} = 'CH';
        elseif  C(x)>0 && H(x)>0 && O(x)==0 && N(x)>0 && S(x)==0 && P(x)==0 && E(x)==0
            t{x,1} = 'CHN';
        elseif  C(x)>0 && H(x)>0 && O(x)==0 && N(x)>0 && S(x)>0 && P(x)==0 && E(x)==0
            t{x,1} = 'CHNS';
        elseif  C(x)>0 && H(x)>0 && O(x)==0 && N(x)>0 && S(x)>0 && P(x)>0 && E(x)==0
            t{x,1} = 'CHNSP';            
        elseif  C(x)>0 && H(x)>0 && O(x)==0 && N(x)>0 && S(x)==0 && P(x)>0 && E(x)==0
            t{x,1} = 'CHNP';                
        elseif  C(x)>0 && H(x)>0 && O(x)==0 && N(x)==0 && S(x)>0 && P(x)>0 && E(x)==0
            t{x,1} = 'CHSP';                
        elseif  C(x)>0 && H(x)>0 && O(x)==0 && N(x)==0 && S(x)>0 && P(x)==0 && E(x)==0
            t{x,1} = 'CHS';
        elseif  C(x)>0 && H(x)>0 && O(x)==0 && N(x)==0 && S(x)==0 && P(x)>0 && E(x)==0
            t{x,1} = 'CHP';    
        elseif  C(x)>0 && H(x)>0 && O(x)>0 && N(x)==0 && S(x)==0 && P(x)==0 && E(x)==0
            t{x,1} = 'CHO';                   
        elseif C(x)>0 && H(x)>0 && O(x)>0 && N(x)>0 && S(x)==0 && P(x)==0 && E(x)==0
            t{x,1} = 'CHON';            
        elseif C(x)>0 && H(x)>0 && O(x)>0 && N(x)==0 && S(x)>0 && P(x)==0 && E(x)==0
            t{x,1} = 'CHOS';           
        elseif C(x)>0 && H(x)>0 && O(x)>0 && N(x)==0 && S(x)==0 && P(x)>0 && E(x)==0
            t{x,1} = 'CHOP';        
        elseif C(x)>0 && H(x)>0 && O(x)>0 && N(x)>0 && S(x)>0 && P(x)==0 && E(x)==0
            t{x,1} = 'CHONS';
        elseif C(x)>0 && H(x)>0 && O(x)>0 && N(x)>0 && S(x)==0 && P(x)>0 && E(x)==0
            t{x,1} = 'CHONP';
        elseif C(x)>0 && H(x)>0 && O(x)>0 && N(x)==0 && S(x)>0 && P(x)>0 && E(x)==0
            t{x,1} = 'CHOSP';
        elseif C(x)>0 && H(x)>0 && O(x)>0 && N(x)>0 && S(x)>0 && P(x)>0 && E(x)==0
            t{x,1} = 'CHONSP';        
        else
            t{x,1} = 'Others';
        end
    end

DataFinal = [num2str(mz,'%4.10f'), num2str(Int), num2str(SN), num2str(EstimC), string(t), num2str(C), num2str(H), num2str(O),...
    num2str(N), num2str(S), num2str(P), num2str(E), num2str(Mass,'%4.10f'),num2str(e),...
    num2str(OC), num2str(HC), num2str(DBE), num2str(DBEC), num2str(AImod),num2str(RelH), string(w)];

DataFinal_Str=string(DataFinal);

%% Exporting Compound Classes to different sheets

Compound_Cl=DataFinal_Str(any(DataFinal_Str=='E',2),:);
Compound_BC=DataFinal_Str(any(DataFinal_Str=='ConAC',2),:);
Compound_SCA=DataFinal_Str(any(DataFinal_Str=='SCA',2),:);
Compound_Lignin=DataFinal_Str(any(DataFinal_Str=='Lignin',2),:);
Compound_Tannin=DataFinal_Str(any(DataFinal_Str=='Tannin',2),:);
Compound_SA=DataFinal_Str(any(DataFinal_Str=='SA',2),:);
Compound_AS=DataFinal_Str(any(DataFinal_Str=='AS',2),:);
Compound_Sugar=DataFinal_Str(any(DataFinal_Str=='Sugar',2),:);
Compound_Prot=DataFinal_Str(any(DataFinal_Str=='Protein',2),:);
Compound_FA=DataFinal_Str(any(DataFinal_Str=='FA',2),:);
Compound_Lip=DataFinal_Str(any(DataFinal_Str=='Lipid',2),:);
Compound_Uns=DataFinal_Str(any(DataFinal_Str=='Unsaturates',2),:);
Compound_Extra=DataFinal_Str(any(DataFinal_Str=='Extra',2),:);

Compound_CHO=DataFinal_Str(any(DataFinal_Str=='CHO',2),:);
Compound_CHON=DataFinal_Str(any(DataFinal_Str=='CHON',2),:);
Compound_CHOP=DataFinal_Str(any(DataFinal_Str=='CHOP',2),:);
Compound_CHOS=DataFinal_Str(any(DataFinal_Str=='CHOS',2),:);
Compound_CHONS=DataFinal_Str(any(DataFinal_Str=='CHONS',2),:);
Compound_CHONP=DataFinal_Str(any(DataFinal_Str=='CHONP',2),:);
Compound_CHOSP=DataFinal_Str(any(DataFinal_Str=='CHOSP',2),:);
Compound_CHONSP=DataFinal_Str(any(DataFinal_Str=='CHONSP',2),:);
Compound_Others=DataFinal_Str(any(DataFinal_Str=='Others',2),:);

%% Statistics 
TotalPeaks=size(Int,1);       
TotalInt=sum(Int);                 
output1=['# of peaks', num2cell(TotalPeaks);'Total Magnitude', num2cell(TotalInt)];

% Statistics for each type of compound class
Compound_Cl_num=str2double(Compound_Cl);
Compound_BC_num=str2double(Compound_BC);
Compound_SCA_num=str2double(Compound_SCA);
Compound_Lignin_num=str2double(Compound_Lignin);
Compound_Tannin_num=str2double(Compound_Tannin);
Compound_SA_num=str2double(Compound_SA);
Compound_AS_num=str2double(Compound_AS);
Compound_Sugar_num=str2double(Compound_Sugar);
Compound_Prot_num=str2double(Compound_Prot);
Compound_FA_num=str2double(Compound_FA);
Compound_Lip_num=str2double(Compound_Lip);
Compound_Uns_num=str2double(Compound_Uns);
Compound_Extra_num=str2double(Compound_Extra);

Compound_CHO_num=str2double(Compound_CHO);
Compound_CHON_num=str2double(Compound_CHON);
Compound_CHOS_num=str2double(Compound_CHOS);
Compound_CHOP_num=str2double(Compound_CHOP);
Compound_CHONS_num=str2double(Compound_CHONS);
Compound_CHONP_num=str2double(Compound_CHONP);
Compound_CHOSP_num=str2double(Compound_CHOSP);
Compound_CHONSP_num=str2double(Compound_CHONSP);
Compound_Others_num=str2double(Compound_Others);
 
[Peaks_Cl,~]=size(Compound_Cl_num);
[Peaks_BC,~]=size(Compound_BC_num);
[Peaks_SCA,~]=size(Compound_SCA_num);
[Peaks_Lignin,~]=size(Compound_Lignin_num);
[Peaks_Tannin,~]=size(Compound_Tannin_num);
[Peaks_SA,~]=size(Compound_SA_num);
[Peaks_AS,~]=size(Compound_AS_num);
[Peaks_Sugar,~]=size(Compound_Sugar_num);
[Peaks_Prot,~]=size(Compound_Prot_num);
[Peaks_FA,~]=size(Compound_FA_num);
[Peaks_Lip,~]=size(Compound_Lip_num);
[Peaks_Uns,~]=size(Compound_Uns_num);
[Peaks_Extra,~]=size(Compound_Extra_num);

[Peaks_CHO,~]=size(Compound_CHO_num);
[Peaks_CHON,~]=size(Compound_CHON_num);
[Peaks_CHOS,~]=size(Compound_CHOS_num);
[Peaks_CHOP,~]=size(Compound_CHOP_num);
[Peaks_CHONS,~]=size(Compound_CHONS_num);
[Peaks_CHONP,~]=size(Compound_CHONP_num);
[Peaks_CHOSP,~]=size(Compound_CHOSP_num);
[Peaks_CHONSP,~]=size(Compound_CHONSP_num);
[Peaks_Others,~]=size(Compound_Others_num);

Int_Cl=sum(Compound_Cl_num(:,2));
Int_BC=sum(Compound_BC_num(:,2));
Int_SCA=sum(Compound_SCA_num(:,2));
Int_Lignin=sum(Compound_Lignin_num(:,2));
Int_Tannin=sum(Compound_Tannin_num(:,2));
Int_SA=sum(Compound_SA_num(:,2));
Int_AS=sum(Compound_AS_num(:,2));
Int_Sugar=sum(Compound_Sugar_num(:,2));
Int_Prot=sum(Compound_Prot_num(:,2));
Int_FA=sum(Compound_FA_num(:,2));
Int_Lip=sum(Compound_Lip_num(:,2));
Int_Uns=sum(Compound_Uns_num(:,2));
Int_Extra=sum(Compound_Extra_num(:,2));

Int_CHO=sum(Compound_CHO_num(:,2));
Int_CHON=sum(Compound_CHON_num(:,2));
Int_CHOS=sum(Compound_CHOS_num(:,2));
Int_CHOP=sum(Compound_CHOP_num(:,2));
Int_CHONS=sum(Compound_CHONS_num(:,2));
Int_CHONP=sum(Compound_CHONP_num(:,2));
Int_CHOSP=sum(Compound_CHOSP_num(:,2));
Int_CHONSP=sum(Compound_CHONSP_num(:,2));
Int_Others=sum(Compound_Others_num(:,2));

% Percentages
PerCHOnum=Peaks_CHO*100/TotalPeaks;
PerCHONnum=Peaks_CHON*100/TotalPeaks;
PerCHOSnum=Peaks_CHOS*100/TotalPeaks;
PerCHOPnum=Peaks_CHOP*100/TotalPeaks;
PerCHONSnum=Peaks_CHONS*100/TotalPeaks;
PerCHONPnum=Peaks_CHONP*100/TotalPeaks;
PerCHOSPnum=Peaks_CHOSP*100/TotalPeaks;
PerCHONSPnum=Peaks_CHONSP*100/TotalPeaks;
PerOthersnum=Peaks_Others*100/TotalPeaks;
PerClnum=Peaks_Cl*100/TotalPeaks;
Sum1=round((PerCHOnum+PerCHONnum+PerCHOSnum+PerCHOPnum+PerCHONSnum+PerCHONPnum+...
    PerCHOSPnum+PerCHONSPnum+PerClnum+PerOthersnum),2);

PerCHOint=Int_CHO*100/TotalInt;
PerCHONint=Int_CHON*100/TotalInt;
PerCHOSint=Int_CHOS*100/TotalInt;
PerCHOPint=Int_CHOP*100/TotalInt;
PerCHONSint=Int_CHONS*100/TotalInt;
PerCHONPint=Int_CHONP*100/TotalInt;
PerCHOSPint=Int_CHOSP*100/TotalInt;
PerCHONSPint=Int_CHONSP*100/TotalInt;
PerOthersint=Int_Others*100/TotalInt;
PerClint=Int_Cl*100/TotalInt;
Sum2=round((PerCHOint+PerCHONint+PerCHOSint+PerCHOPint+PerCHONSint+PerCHONPint+...
    PerCHOSPint+PerCHONSPint+PerClint+PerOthersint),2);

PerBCnum=Peaks_BC*100/TotalPeaks;
PerSCAnum=Peaks_SCA*100/TotalPeaks;
PerLigninnum=Peaks_Lignin*100/TotalPeaks;
PerTanninnum=Peaks_Tannin*100/TotalPeaks;
PerSAnum=Peaks_SA*100/TotalPeaks;
PerASnum=Peaks_AS*100/TotalPeaks;
PerSugarnum=Peaks_Sugar*100/TotalPeaks;
PerProtnum=Peaks_Prot*100/TotalPeaks;
PerFAnum=Peaks_FA*100/TotalPeaks;
PerLipnum=Peaks_Lip*100/TotalPeaks;
PerUnsnum=Peaks_Uns*100/TotalPeaks;
PerExtranum=Peaks_Extra*100/TotalPeaks;

Sum3=round((PerBCnum+PerSCAnum+PerLigninnum+PerTanninnum+PerSAnum+PerASnum+...
    PerSugarnum+PerProtnum+PerFAnum+PerLipnum+PerUnsnum+PerExtranum+PerClnum),2);

PerBCint=Int_BC*100/TotalInt;
PerSCAint=Int_SCA*100/TotalInt;
PerLigninint=Int_Lignin*100/TotalInt;
PerTanninint=Int_Tannin*100/TotalInt;
PerSAint=Int_SA*100/TotalInt;
PerASint=Int_AS*100/TotalInt;
PerSugarint=Int_Sugar*100/TotalInt;
PerProtint=Int_Prot*100/TotalInt;
PerFAint=Int_FA*100/TotalInt;
PerLipint=Int_Lip*100/TotalInt;
PerUnsint=Int_Uns*100/TotalInt;
PerExtraint=Int_Extra*100/TotalInt;

Sum4=round((PerBCint+PerSCAint+PerLigninint+PerTanninint+PerSAint+PerASint+...
    PerSugarint+PerProtint+PerFAint+PerLipint+PerUnsint+PerExtraint+PerClint),2);

if Sum1~=100
    error('Error! The sum1 of percentages does not equal 100%!')
end

if Sum2~=100
    error('Error! The sum2 of percentages does not equal 100%!')
end

if Sum3~=100
    error('Error! The sum3 of percentages does not equal 100%!')
end

if Sum4~=100
    error('Error! The sum4 of percentages does not equal 100%!')
end

%% Figures

figurename=char(filename(1:end-5));
figure(1) 
set(gcf,'WindowState','maximized')
    subplot(2,2,1) % vK all
        hold on
        title('vK All Formulas');
        xlabel('O/C'); ylabel('H/C');
        xlim([0 1.2]); ylim([0 2.5]);
        plot(OC,HC,'r.','MarkerSize',5);
        h1=refline(-1,2);           % AImod=0
        h2=refline(-0.3845,1.0584); % AImod=0.5
        h3=refline(-0.3740,0.7551); % AImod=0.67
        h1.Color='k';h2.Color='k';h3.Color='k';
        h1.LineWidth=1.5;h2.LineWidth=1.5;h3.LineWidth=1.5;
        set(gcf,'color',[0.85 0.85 0.85]);
        text(0.2,2.4,figurename,'FontSize',15,'FontWeight','bold','Interpreter','none')
        axis manual
        set(gca,'fontweight','bold') 
        set(gca,'FontSize',15)
        hold off
    subplot(2,2,[2,4]) 
        title('3D vK with Magnitude');
        hold on
        IntLog=log(RelH);
        scatter3(OC(:),HC(:),IntLog(:),70,IntLog(:),'filled','o');
        colormap(jet);
        handle_colorbar=colorbar;
        set(get(handle_colorbar,'title'),'string','log(Rel.Magn.)','fontsize',13);
        set(handle_colorbar,'Position',[0.939895833333331,0.108961578400831,0.011111111111111,0.815]);
        view(45,30)
        grid on        
        xlabel('O/C') % x-axid label
        ylabel('H/C') % y-axis label
        zlabel('Rel. Abundance (in log units)')
        xlim([0 1.2]); ylim([0 2.5]);
        axis manual
        set(gca,'fontweight','bold')
        set(gca,'FontSize',15)
        hold off
    subplot(2,2,3)
        hold on
        title('H/C vs Molecular Weight');
        xlabel('Molecular Weight'); ylabel('H/C');
        xlim([0 1.2]); ylim([0 2.5]);
        plot(Mass,HC,'r.','MarkerSize',5);
        axis manual
        xlim([150 1000]); ylim([0 2.5]);
        xticks([200 400 600 800 1000]);
        set(gca,'fontweight','bold')
        set(gca,'FontSize',15)
        hold off
print(gcf,['vK_all_' figurename],'-dpng','-r300');

figure(2)
     subplot(121) % vK different types (eg. CHO, CHON,...)
        hold on
        set(gcf,'color',[0.85 0.85 0.85]);
        title('vK Formula Types');
        xlabel('O/C'); ylabel('H/C');
        xlim([0 1.2]); ylim([0 2.5]);
        OC_CHO=Compound_CHO_num(:,format.Column_OC); HC_CHO=Compound_CHO_num(:,format.Column_HC);
        OC_CHON=Compound_CHON_num(:,format.Column_OC); HC_CHON=Compound_CHON_num(:,format.Column_HC);
        OC_CHOS=Compound_CHOS_num(:,format.Column_OC); HC_CHOS=Compound_CHOS_num(:,format.Column_HC);
        OC_CHOP=Compound_CHOP_num(:,format.Column_OC); HC_CHOP=Compound_CHOP_num(:,format.Column_HC);
        OC_CHONS=Compound_CHONS_num(:,format.Column_OC); HC_CHONS=Compound_CHONS_num(:,format.Column_HC);
        OC_CHONP=Compound_CHONP_num(:,format.Column_OC); HC_CHONP=Compound_CHONP_num(:,format.Column_HC);
        OC_CHOSP=Compound_CHOSP_num(:,format.Column_OC); HC_CHOSP=Compound_CHOSP_num(:,format.Column_HC);
        OC_CHONSP=Compound_CHONSP_num(:,format.Column_OC); HC_CHONSP=Compound_CHONSP_num(:,format.Column_HC);
        OC_Others=Compound_Others_num(:,format.Column_OC); HC_Others=Compound_Others_num(:,format.Column_HC);
        OC_Cl=Compound_Cl_num(:,format.Column_OC); HC_Cl=Compound_Cl_num(:,format.Column_HC);
        plot(OC_CHO,HC_CHO,'r.','MarkerSize',13,'DisplayName', 'CHO');
        plot(OC_CHON,HC_CHON,'b.','MarkerSize',13,'DisplayName', 'CHON');
        plot(OC_CHOS,HC_CHOS,'k.','MarkerSize',13,'DisplayName', 'CHOS');
        plot(OC_CHOP,HC_CHOP,'m.','MarkerSize',13,'DisplayName', 'CHOP');
        plot(OC_CHONS,HC_CHONS,'b*','MarkerSize',13,'DisplayName', 'CHONS');
        plot(OC_CHONP,HC_CHONP,'b+','MarkerSize',13,'DisplayName', 'CHONP');
        plot(OC_CHOSP,HC_CHOSP,'k*','MarkerSize',13,'DisplayName', 'CHOSP');
        plot(OC_CHONSP,HC_CHONSP,'bh','MarkerSize',13,'DisplayName', 'CHONSP');
        plot(OC_Cl,HC_Cl,'g.','MarkerSize',13,'DisplayName', 'E');
        plot(OC_Others,HC_Others,'bh','MarkerSize',13,'DisplayName', 'Others');
        set(gca,'fontweight','bold')
        set(gca,'ycolor','k')
        set(gca,'xcolor','k')
        legend show
        h1=refline(-1,2);           % AImod=0    
        h2=refline(-0.3845,1.0584); % AImod=0.5
        h3=refline(-0.3740,0.7551); % AImod=0.67
        h1.Color='k';h2.Color='k';h3.Color='k';
        h1.LineWidth=1.5;h2.LineWidth=1.5;h3.LineWidth=1.5;
        h1.DisplayName='off';h2.DisplayName='off';h3.DisplayName='off';
        h1.HandleVisibility='off';h2.HandleVisibility='off';h3.HandleVisibility='off';
        l1=legend('Location','northeast'); 
        set(l1,'NumColumns',2);
        text(0.2,2.4,figurename,'FontSize',15,'FontWeight','bold','Interpreter','none')
        axis manual
        set(gca,'FontSize',15)
        set(l1,'FontSize',9)
        hold off 
     subplot(122) % vK different compound classes, e.g., lignin, BC, etc. 
        hold on
        title('vK Compound Classes');
        xlabel('O/C'); ylabel('H/C');
        xlim([0 1.2]); ylim([0 2.5]);
        OC_BC=Compound_BC_num(:,format.Column_OC); HC_BC=Compound_BC_num(:,format.Column_HC);
        OC_SCA=Compound_SCA_num(:,format.Column_OC); HC_SCA=Compound_SCA_num(:,format.Column_HC);
        OC_Lignin=Compound_Lignin_num(:,format.Column_OC); HC_Lignin=Compound_Lignin_num(:,format.Column_HC);
        OC_Tannin=Compound_Tannin_num(:,format.Column_OC); HC_Tannin=Compound_Tannin_num(:,format.Column_HC);
        OC_SA=Compound_SA_num(:,format.Column_OC); HC_SA=Compound_SA_num(:,format.Column_HC);
        OC_AS=Compound_AS_num(:,format.Column_OC); HC_AS=Compound_AS_num(:,format.Column_HC);
        OC_Sugar=Compound_Sugar_num(:,format.Column_OC); HC_Sugar=Compound_Sugar_num(:,format.Column_HC);
        OC_Prot=Compound_Prot_num(:,format.Column_OC); HC_Prot=Compound_Prot_num(:,format.Column_HC);
        OC_FA=Compound_FA_num(:,format.Column_OC); HC_FA=Compound_FA_num(:,format.Column_HC);
        OC_Lip=Compound_Lip_num(:,format.Column_OC); HC_Lip=Compound_Lip_num(:,format.Column_HC);
        OC_Uns=Compound_Uns_num(:,format.Column_OC); HC_Uns=Compound_Uns_num(:,format.Column_HC);
        OC_Extra=Compound_Extra_num(:,format.Column_OC); HC_Extra=Compound_Extra_num(:,format.Column_HC);
        plot(OC_BC,HC_BC,'Marker','*','Color','k','MarkerSize',13,'LineStyle','none','DisplayName','ConAC');
        plot(OC_Lignin,HC_Lignin,'Marker','.','Color',[0.07,0.62,1.00],'MarkerSize',13,'LineStyle','none','DisplayName','Lignin');
        plot(OC_Tannin,HC_Tannin,'Marker','.','Color',[0,0,1.00],'MarkerSize',13,'LineStyle','none','DisplayName','Tannin');
        plot(OC_SCA,HC_SCA,'Marker','.','Color',[1.00,0.07,0.65],'MarkerSize',24,'LineStyle','none','DisplayName','SCA');
        plot(OC_SA,HC_SA,'Marker','x','Color',[0,0,0],'MarkerSize',24,'LineStyle','none','DisplayName','Sulfonic Acid');
        plot(OC_AS,HC_AS,'Marker','.','Color',[0,1,1],'MarkerSize',24,'LineStyle','none','DisplayName','Amino Sugar');
        plot(OC_Sugar,HC_Sugar,'k.','MarkerSize',24,'LineStyle','none','DisplayName','Sugar');
        plot(OC_Prot,HC_Prot,'Marker','.','Color',[0.00,1.00,0.07],'MarkerSize',24,'LineStyle','none','DisplayName','Protein');
        plot(OC_FA,HC_FA,'Marker','pentagram','Color',[1.00,0.41,0.16],'MarkerSize',13,'LineStyle','none','DisplayName','Fatty Acid');
        plot(OC_Lip,HC_Lip,'.r','MarkerSize',24,'LineStyle','none','DisplayName','Lipid');
        plot(OC_Uns,HC_Uns,'Marker','pentagram','Color',[0.06,1.00,1.00],'MarkerSize',24,'LineStyle','none','DisplayName','Unsaturates');
        plot(OC_Cl,HC_Cl,'Marker','pentagram','Color',[0.72,0.27,1.00],'MarkerSize',24,'LineStyle','none','DisplayName','E');
        plot(OC_Extra,HC_Extra,'kh','MarkerSize',13,'LineStyle','none','DisplayName','Extra');
        h1=refline(-1,2);           % AImod=0
        h2=refline(-0.3845,1.0584); % AImod=0.5
        h3=refline(-0.3740,0.7551); % AImod=0.67
        h1.Color='k';h2.Color='k';h3.Color='k';
        h1.LineWidth=1.5;h2.LineWidth=1.5;h3.LineWidth=1.5;
        h1.HandleVisibility='off';h2.HandleVisibility='off';h3.HandleVisibility='off';
        l2=legend('Location','south'); 
        set(l2,'NumColumns',4);
        text(0.2,2.4,figurename,'FontSize',15,'FontWeight','bold','Interpreter','none')
        set(gca,'fontweight','bold')
        set(gca,'ycolor','k')
        set(gca,'xcolor','k')
        axis manual
        set(gca,'FontSize',15)
        set(l2,'FontSize',9)
        hold off
    set(gcf,'WindowState','maximized')
    print(gcf,['vK_groups_' figurename],'-dpng','-r300');

%% Export

filename_final=[filename(1:end-5) '_CompoundClass.xlsx'];
warning('off','MATLAB:xlswrite:AddSheet')

xlswrite(filename_final,[titles,'RelMagn','Class';DataFinal],'Sheet1','A1'); % Exports original data

if PerCHOnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_CHO],'CHO','A1');
end

if PerCHONnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_CHON],'CHON','A1');
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_CHO;Compound_CHON],'CHO+CHON','A1')
end

if PerCHOPnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_CHOP],'CHOP','A1');
end

if PerCHOSnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_CHOS],'CHOS','A1');
end

if PerCHONSnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_CHONS],'CHONS','A1');
end

if PerCHONPnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_CHONP],'CHONP','A1');
end

if PerCHOSPnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_CHOSP],'CHOSP','A1');
end

if PerCHONSPnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_CHONSP],'CHONSP','A1');
end

if PerOthersnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_Others],'Other types','A1');
end

if PerClnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_Cl],'E','A1');
end

if PerBCnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_BC],'ConAC','A1');
end

if PerSCAnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_SCA],'SCA','A1');
end

if PerLigninnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_Lignin],'Lignin','A1');
end

if PerTanninnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_Tannin],'Tannin','A1');
end

if PerSAnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_SA],'Sulfonic Acid','A1');
end

if PerASnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_AS],'Amino Sugar','A1');
end

if PerSugarnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_Sugar],'Sugar','A1');
end

if PerProtnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_Prot],'Protein','A1');
end

if PerFAnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_FA],'Fatty Acid','A1');
end

if PerLipnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_Lip],'Lipids','A1');
end

if PerUnsnum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_Uns],'Unsaturates','A1');
end

if PerExtranum >0
    xlswrite(filename_final,[titles,'RelMagn','Class';Compound_Extra],'Extra','A1');
end

titles_percentages1={'CHO'; 'CHON';'CHO+CHON'; 'CHOS'; 'CHOP'; 'CHONS'; 'CHONP'; 'CHOSP'; 'CHONSP'; 'E'; 'Others'};
titles_percentages2={'ConAC';'SCA';'Lignin';'Tannin';'Sulfonic Acid';'Amino Sugar';'Sugar';...
    'Protein';'Fatty Acid';'Lipid';'Unsaturates';'E';'Extra';};
titles_percentages3={'# formulas','Num%','Magn%'};

peaks_CHOCHON=Peaks_CHO+Peaks_CHON;
PerCHOCHONnum=PerCHOnum+PerCHONnum;
PerCHOCHONint=PerCHOint+PerCHONint;

stats1=[Peaks_CHO,round(PerCHOnum,2),round(PerCHOint,2);Peaks_CHON,...
    round(PerCHONnum,2),round(PerCHONint,2);peaks_CHOCHON,...
    round(PerCHOCHONnum,2),round(PerCHOCHONint,2);Peaks_CHOS,round(PerCHOSnum,2),...
    round(PerCHOSint,2);Peaks_CHOP,round(PerCHOPnum,2),round(PerCHOPint,2);...
    Peaks_CHONS,round(PerCHONSnum,2),round(PerCHONSint,2);Peaks_CHONP,...
    round(PerCHONPnum,2),round(PerCHONPint,2);Peaks_CHOSP,round(PerCHOSPnum,2),...
    round(PerCHOSPint,2);Peaks_CHONSP,round(PerCHONSPnum,2),round(PerCHONSPint,2);...
    Peaks_Cl,round(PerClnum,2),round(PerClint,2);Peaks_Others,round(PerOthersnum,2),round(PerOthersint,2)];

stats2=[Peaks_BC,round(PerBCnum,2),round(PerBCint,2);Peaks_SCA,round(PerSCAnum,2),...
    round(PerSCAint,2);Peaks_Lignin,round(PerLigninnum,2),round(PerLigninint,2);Peaks_Tannin,round(PerTanninnum,2),...
    round(PerTanninint,2);Peaks_SA,round(PerSAnum,2),round(PerSAint,2);Peaks_AS,...
    round(PerASnum,2),round(PerASint,2);Peaks_Sugar,round(PerSugarnum,2),round(PerSugarint,2);...
    Peaks_Prot,round(PerProtnum,2),round(PerProtint,2);Peaks_FA,round(PerFAnum,2),round(PerFAint,2);...
    Peaks_Lip,round(PerLipnum,2),round(PerLipint,2);Peaks_Uns,round(PerUnsnum,2),round(PerUnsint,2);...   
    Peaks_Cl,round(PerClnum,2),round(PerClint,2);Peaks_Extra,round(PerExtranum,2),round(PerExtraint,2)];
 
output2=[titles_percentages1,string(stats1)];
output3=[titles_percentages2,string(stats2)];

xlswrite(filename_final,output1,'Sheet1','W2');
xlswrite(filename_final,[' ', titles_percentages3;output2],'Sheet1','W5');
xlswrite(filename_final,[' ', titles_percentages3;output3],'Sheet1','W18');


%% Perform FTMS_Metrics on all

FTMS_Metrics(filename_final,'Sheet1','W33','W37',true);
loc1='W3'; loc2='W7';

if PerCHOnum >0
    FTMS_Metrics(filename_final,'CHO',loc1,loc2,true)
end

if PerCHONnum >0
    FTMS_Metrics(filename_final,'CHON',loc1,loc2,true)
    FTMS_Metrics(filename_final,'CHO+CHON',loc1,loc2,true)
end

if PerCHOPnum >0
    FTMS_Metrics(filename_final,'CHOP',loc1,loc2,true)
end

if PerCHOSnum >0
    FTMS_Metrics(filename_final,'CHOS',loc1,loc2,true)
end

if PerCHONSnum >0
    FTMS_Metrics(filename_final,'CHONS',loc1,loc2,true)
end

if PerCHONPnum >0
    FTMS_Metrics(filename_final,'CHONP',loc1,loc2,true)
end

if PerCHOSPnum >0
    FTMS_Metrics(filename_final,'CHOSP',loc1,loc2,true)
end

if PerCHONSPnum >0
    FTMS_Metrics(filename_final,'CHONSP',loc1,loc2,true)
end

if PerOthersnum >0
    FTMS_Metrics(filename_final,'Other types',loc1,loc2,true)
end

if PerClnum >0
    FTMS_Metrics(filename_final,'E',loc1,loc2,true)
end

if PerBCnum >0
    FTMS_Metrics(filename_final,'ConAC',loc1,loc2,true)
end

if PerSCAnum >0
    FTMS_Metrics(filename_final,'SCA',loc1,loc2,true)
end

if PerLigninnum >0
    FTMS_Metrics(filename_final,'Lignin',loc1,loc2,true)
end

if PerTanninnum >0
    FTMS_Metrics(filename_final,'Tannin',loc1,loc2,true)
end

if PerSAnum >0
    FTMS_Metrics(filename_final,'Sulfonic Acid',loc1,loc2,true)
end

if PerASnum >0
    FTMS_Metrics(filename_final,'Amino Sugar',loc1,loc2,true)
end

if PerSugarnum >0
    FTMS_Metrics(filename_final,'Sugar',loc1,loc2,true)
end

if PerProtnum >0
    FTMS_Metrics(filename_final,'Protein',loc1,loc2,true)
end

if PerFAnum >0
    FTMS_Metrics(filename_final,'Fatty Acid',loc1,loc2,true)
end

if PerLipnum >0
    FTMS_Metrics(filename_final,'Lipids',loc1,loc2,true)
end

if PerUnsnum >0
    FTMS_Metrics(filename_final,'Unsaturates',loc1,loc2,true)
end

if PerExtranum >0
    FTMS_Metrics(filename_final,'Extra',loc1,loc2,true)
end

Stats=[TotalPeaks;stats1(:,1);stats1(:,2);stats1(:,3);stats2(:,1);stats2(:,2);stats2(:,3)];
if nargout == 1
    varargout{1}=[filename(1:end-5);string(Stats)];
else
    ;
end

disp(['Finished categorizing the compound classes of ' char(filename)])
end