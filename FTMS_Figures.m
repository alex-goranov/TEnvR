function FTMS_Figures(filename)
%% Description

% This code loads a list of formulas and plots a variety of figures

% Example:
% FTMS_Figures('Sample1_Final.xlsx')

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
format = FTMS_ConfigurationToolbox;

[Data,TXT]=xlsread(filename);
titles=TXT(1,:);

Data=sortrows(Data,format.Column_ExactMass);
mz=Data(:,format.Column_mz);
Int=Data(:,format.Column_Magnitude);
C=Data(:,format.Column_C);
H=Data(:,format.Column_H);
O=Data(:,format.Column_O);
N=Data(:,format.Column_N);
S=Data(:,format.Column_S);
P=Data(:,format.Column_P);
E=Data(:,format.Column_E);

e=Data(:,format.Column_error);
OC=Data(:,format.Column_OC);
HC=Data(:,format.Column_HC);
DBE=Data(:,format.Column_DBE);
DBEC=Data(:,format.Column_DBEC);
AImod=Data(:,format.Column_AImod);

Data(:,format.Column_ExactMass)=Data(:,format.Column_C)*format.Mass_12C+Data(:,format.Column_H)*format.Mass_1H+Data(:,format.Column_N)*format.Mass_14N+...
    Data(:,format.Column_O)*format.Mass_16O+Data(:,format.Column_S)*format.Mass_32S+Data(:,format.Column_P)*format.Mass_31P+Data(:,format.Column_E)*format.Mass_Heteroelement;
Data(:,format.Column_ExactMass)=round(Data(:,format.Column_ExactMass),format.Precision+1); % Round to the uncertainty in the original data

Mass=Data(:,format.Column_ExactMass);

global Precision
Precision=format.Precision;

NC=N./C;   
SC=S./C;
PC=P./C;
ClC=E./C;
DBEH=DBE./H;
DBEO=DBE./O;
DBEminusO=DBE-O;

if format.Heteroelement_Halogen
    NOSC=4-(((4*C)+H-(3*N)-(2*O)-(2*S)+(5*P)-E)./C);
elseif format.Heteroelement_N15
	NOSC=4-(((4*C)+H-(3*N)-(3*E)-(2*O)-(2*S)+(5*P))./C);
else
    NOSC=4-(((4*C)+H-(3*N)-(2*O)-(2*S)+(5*P))./C);
end

% Finding condensed aromatic molecules "black carbon" (AImod >= 0.67, C >= 15)
BlackCarbon=Data(AImod >= 0.67 & C >= 15,format.Column_mz:format.Column_AImod);
BlackCarbon_DBE=BlackCarbon(:,format.Column_DBE);
BlackCarbon_C=BlackCarbon(:,format.Column_C);
BlackCarbon_Int=BlackCarbon(:,format.Column_Magnitude);
Ring=(BlackCarbon_DBE-3.1374)/2.436;
RingC=Ring./BlackCarbon_C;

IntRel=Int/sum(Int);
IntLog=log(IntRel);
BlackCarbon_IntRel=BlackCarbon_Int/sum(Int);
BlackCarbon_IntLog=log(BlackCarbon_IntRel);

figure % Histograms
set(gcf,'WindowState','maximized');
    subplot(3,8,1) 
        his(C,'C-number Histogram','C-number');
    subplot(3,8,2)
        his(H,'H-number Histogram','H-number');
    subplot(3,8,3)
        his(O,'O-number Histogram','O-number');
    subplot(3,8,4)
        his(N,'N-number Histogram','N-number');
    subplot(3,8,5)
        his(S,'S-number Histogram','S-number');
    subplot(3,8,6)
        his(P,'P-number Histogram','P-number');
    subplot(3,8,7)
        his(E,'E-number Histogram','E-number');
    subplot(3,8,8)
        his(Mass,'ExactMass Histogram','ExactMass (Da)');
    subplot(3,8,9) 
        his(e,'Error Histogram','Error (ppm)');
    subplot(3,8,10)
        his(HC,'H/C Histogram','H/C Ratio');
    subplot(3,8,11)
        his(OC,'O/C Histogram','O/C Ratio');
    subplot(3,8,12)
        his(NC,'N/C Histogram','N/C Ratio');
    subplot(3,8,13)
        his(SC,'S/C Histogram','S/C Ratio');
    subplot(3,8,14)
        his(PC,'P/C Histogram','P/C Ratio');
    subplot(3,8,15)
        his(ClC,'E/C Histogram','E/C Ratio');
    subplot(3,8,16)
        his(NOSC,'NOSC Histogram','Nominal Oxidation State of Carbon');
    subplot(3,8,17)
        his(DBE,'DBE Histogram','Double Bond Equivalency (DBE)');
    subplot(3,8,18)
        his(DBEC,'DBE/C Histogram','C-Normalized DBE');
    subplot(3,8,19)
        his(DBEH,'DBE/H Histogram','H-Normalized DBE');
    subplot(3,8,20)
        his(DBEO,'DBE/O Histogram','O-Normalized DBE');
    subplot(3,8,21)
        his(DBEminusO,'DBE-O Histogram','O-Corrected DBE');
    subplot(3,8,22)
        his(AImod,'AI_{MOD} Histogram','Modified Aromaticity Index (AI_{MOD})');
    subplot(3,8,23)
        his(Ring,'Ring Histogram','# of Aromatic Rings');
    subplot(3,8,24)
        his(RingC,'Ring/C Histogram','C-Normalized # of Aromatic Rings');
print(gcf,['Fig_Histograms_' char(filename(1:end-5)) '.png'],'-dpng','-r300');
    
figure % Flat 3D vKs
set(gcf,'WindowState','maximized');
    subplot(2,4,1)
        ThreeDvK(IntLog,'3D vK w/ Magnitude','Rel. Magnitude (in log units)',OC,HC)
    subplot(2,4,2)
        ThreeDvK(Mass,'3D vK w/ Molecular Weight','Molecular Weight',OC,HC)
    subplot(2,4,3)
        ThreeDvK(NC,'3D vK w/ N/C','N/C Ratio',OC,HC)
    subplot(2,4,4)
        ThreeDvK(SC,'3D vK w/ S/C','S/C Ratio',OC,HC)
    subplot(2,4,5)
        ThreeDvK(PC,'3D vK w/ P/C','P/C Ratio',OC,HC)
    subplot(2,4,6)
        ThreeDvK(ClC,'3D vK w/ E/C','E/C Ratio',OC,HC)
    subplot(2,4,7)
        ThreeDvK(DBE,'3D vK w/ DBE','Double Bond Equivalency',OC,HC)
    subplot(2,4,8)
        ThreeDvK(NOSC,'3D vK w/ NOSC','Nominal Oxidation State of Carbon',OC,HC)
print(gcf,['Fig_3DvKs_' char(filename(1:end-5)) '.png'],'-dpng','-r300');

figure % KMD plots
set(gcf,'WindowState','maximized');
    subplot(2,4,1)
        KMD(Mass,14.01565,14,'CH_2') 
    subplot(2,4,2)    
        KMD(Mass,2.01565,2,'H_2') 
    subplot(2,4,3)    
        KMD(Mass,43.98983,44,'COO') 
    subplot(2,4,4)    
        KMD(Mass,44.026216,44,'CH_2O') 
    subplot(2,4,5)    
        KMD(Mass,31.9898292,32,'O_2') 
    subplot(2,4,6)    
        KMD(Mass,18.0105646,18,'H_2O') 
    subplot(2,4,7)    
        KMD(Mass,17.026549,17,'NH_3') 
    subplot(2,4,8)    
        KMD(Mass,35.9766777,36,'HCl') 
print(gcf,['Fig_KMDplots_' char(filename(1:end-5)) '.png'],'-dpng','-r300');

figure % C-number plots
set(gcf,'WindowState','maximized');        
    subplot(2,5,1) 
        Cnumber(C,HC,IntLog,'H/C Ratio')
    subplot(2,5,2)
        Cnumber(C,OC,IntLog,'O/C Ratio')
    subplot(2,5,3)
        Cnumber(C,NC,IntLog,'N/C Ratio')
    subplot(2,5,4)
        Cnumber(C,SC,IntLog,'S/C Ratio')
    subplot(2,5,5)
        Cnumber(C,PC,IntLog,'P/C Ratio')
    subplot(2,5,6)
        Cnumber(C,NOSC,IntLog,'Nominal Oxidation State of Carbon')
    subplot(2,5,7)
        Cnumber(C,DBE,IntLog,'Double Bond Equivalency (DBE)')
    subplot(2,5,8)
        Cnumber(C,DBEC,IntLog,'C-Normalized DBE')
    subplot(2,5,9) 
        Cnumber(C,AImod,IntLog,'Modified Aromaticity Index (AI_{MOD})')
    subplot(2,5,10)
        Cnumber(BlackCarbon_C,RingC,BlackCarbon_IntLog,'C-Normalized # of Aromatic Rings')
        
print(gcf,['Fig_Cnumberplots_' char(filename(1:end-5)) '.png'],'-dpng','-r300');

disp(['Finished plotting figures of ' char(filename)])
end

%% Internal Functions

% Plotting a histogram with 10 bins
function his(X,titl,xlabl)
    hold on
    h=histogram(X,10); % 10 bins
    title(titl);
    xlabel(xlabl); ylabel('Number of Formulas');
    x=gca;
    axis manual 
    set(h,'FaceColor',[0.00,0.45,0.74])
    x.FontWeight='bold'; 
    hold off
end
    
% Plotting a 3D vK
function ThreeDvK(var,titl,zlabl,OC,HC)
    hold on
    title(titl);
    scatter3(OC(:),HC(:),var(:),50,var(:),'filled','o');
    colormap(jet);
    colorbar;
    %view(-27,49)
    grid on        
    xlabel('O/C') 
    ylabel('H/C') 
    zlabel(zlabl)
    xlim([0 1.2]); ylim([0 2.5]);
    x=gca;
    x.FontWeight='bold'; 
    axis manual
    hold off    
end

% KMD series plot
function KMD(Mass,SeriesMass,SeriesNominalMass,label) 
    global Precision
    KM=round(Mass*(SeriesNominalMass/SeriesMass));
    KMD=Mass*(SeriesNominalMass/SeriesMass)-KM;
    KMD_round=round(KMD,Precision);
    hold on
    t=scatter(KM,KMD_round,200,'.k');
    xlabel(['Kendrick Mass (' label ')'])        
    ylabel(['Kendrick Mass Defect (' label ')']) 
    set(t,'SizeData',70);
    xlim([200 800]); ylim([-1 1]);
    x=gca;
    x.FontWeight='bold'; 
    axis manual
    hold off
end

% C-number plot
function Cnumber(C,variable,IntLog,ylabl)
    hold on
    scatter(C,variable,50,IntLog,'filled','o')
    colormap(jet)
    x=gca;
    x.FontWeight='bold'; 
    xlabel('# of Carbons')  
    ylabel(ylabl)           
    colorbar;
    hold off
end