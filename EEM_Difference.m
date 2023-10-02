function EEM_Difference(filename1,filename2)
%% Description

% This code takes two processed EEM spectra and calculates a difference spectrum.
% The spectra must be located in .csv files (comma-separated values). 

% Examples:
% EEM_Difference('Sample 1_ref_Final.csv','Sample 2_ref_Final.csv')

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

clearvars -except filename1 filename2 
close all

global Peak_Labels % Do not edit this line!

% Peak Labels
Peak_Labels=true; % choose true or false

% EEM orientation: choose if you want the excitation to be as x-axis
Excitation_Xaxis=false; % choose true if you want excitation to be as x-exis; false = y-axis

% Normalization: choose  if you want to normalize the EEMs to their total spectral intensity
Normalization=true; % choose true or false

%% Stage1: Import and organize data

Data1=csvread(filename1); % Import data
Data2=csvread(filename2); % Import data

Ex1=Data1(1,:); Em1=Data1(:,1); Ex2=Data2(1,:); Em2=Data2(:,1); % Extract excitation and emission ranges
Ex1=Ex1(2:end); Em1=Em1(2:end); Ex2=Ex2(2:end); Em2=Em2(2:end); % Extract excitation and emission ranges
Int1=Data1; Int2=Data2; Int1(1,:)=[]; Int1(:,1)=[]; Int2(1,:)=[]; Int2(:,1)=[]; % Extract intensity data 
Int1=transpose(Int1); Int2=transpose(Int2);

% Make sure the dimensions are consistent        
Size_Em1=size(Em1);
Size_Em2=size(Em2);
Size_Ex1=size(Ex1);
Size_Ex2=size(Ex2);
Size_Int1=size(Int1);
Size_Int2=size(Int2);

if Size_Em1(1)~= Size_Int1(2)
    error('Error! Check the dimensions of the emission of sample 1')
end
if Size_Em2(1)~= Size_Int2(2)
    error('Error! Check the dimensions of the emission of sample 2')
end
if Size_Ex1(2)~= Size_Int1(1)
    error('Error! Check the dimensions of the excitation of sample 1')
end
if Size_Ex2(2)~= Size_Int2(1)
    error('Error! Check the dimensions of the excitation of sample 2')
end

%% Stage 2: Compare samples and interpolate if necessary
if Size_Int1(1)~= Size_Int2(1) || Size_Int1(2)~= Size_Int2(2) 
        
    % Creates X and Y matrices of the old EX and EM wavelengths
    [Interp_X1,Interp_Y1] = meshgrid(Ex1,Em1);
    [Interp_X2,Interp_Y2] = meshgrid(Ex2,Em2);
    
    % V is the old Int values
    Interp_V1=Int1;
    Interp_V2=Int2;
    
    % Determine new Ex and EM dimensions
    Ex_Interp=round(min([Ex1,Ex2])):1:round(max([Ex1,Ex2]));
    Em_Interp=round(min([Em1;Em2])):1:round(max([Em1;Em2]));
    Em_Interp=transpose(Em_Interp);

    % Creates new X and Y matrices of the old EX and EM wavelengths
    [Interp_Xq,Interp_Yq] = meshgrid(Ex_Interp,Em_Interp);

    % Interpoloate
    Int_Interp1=interp2(Interp_X1,Interp_Y1,Interp_V1',Interp_Xq,Interp_Yq);
    Int_Interp1=fillmissing(Int_Interp1,'constant',0); %Any NaN values are turned into 0
    
    Int_Interp2=interp2(Interp_X2,Interp_Y2,Interp_V2',Interp_Xq,Interp_Yq);
    Int_Interp2=fillmissing(Int_Interp2,'constant',0); %Any NaN values are turned into 0

    Data_Interp1=[0,Ex_Interp;Em_Interp,Int_Interp1];
    Data_Interp2=[0,Ex_Interp;Em_Interp,Int_Interp2];
    
    % Change initial variables
    Ex1=Ex_Interp;
    Ex2=Ex_Interp;
    Em1=Em_Interp;
    Em2=Em_Interp;
    Int1=Int_Interp1;
    Int2=Int_Interp2;
end

if Normalization
    Int1=Int1/nansum(nansum(Int1));
    Int2=Int2/nansum(nansum(Int2));
end

% Plotting
figure
    subplot(1,3,1);
        hold on
        if Excitation_Xaxis 
            plotEEM_SwitchedAxis(char(filename1(1:end-4)),Int1,Em1,Ex1)
        else
            plotEEM(char(filename1(1:end-4)),Int1,Em1,Ex1)
        end
        hold off
    subplot(1,3,2);
        hold on
        if Excitation_Xaxis 
            plotEEM_SwitchedAxis(char(filename2(1:end-4)),Int2,Em2,Ex2)
        else
            plotEEM(char(filename2(1:end-4)),Int2,Em2,Ex2)
        end
        hold off
    subplot(1,3,3);
        Diff=Int2-Int1;
        hold on
        if Excitation_Xaxis 
            plotEEM_SwitchedAxis('Difference',Diff,Em1,Ex1)
        else
            plotEEM('Difference',Diff,Em1,Ex1)
        end
        hold off
set(gcf,'Units','Normalized');
set(gcf,'OuterPosition',[-0.003645833333333,0.465740740740741,1.003645833333333,0.517592592592593]);

figurename=[char(filename1(1:end-4)) '_' char(filename2(1:end-4))];
print(gcf,['DifferenceEEM_' figurename '.tiff'],'-dtiff','-r500');

disp(['Subtracted ' char(filename2) ' from ' char(filename1) ' to produce a difference EEM spectrum!'])
end                  

%% Internal Functions

% Plotting EEM with x-axis = emission
function plotEEM(tit,Signal,Em,Ex)
global  Peak_Labels
        title(tit,'Interpreter','none','fontsize',13,'Position',[median(Em),max(Ex)+10,0]); 
        [~,h]=contourf(Em,Ex,squeeze(Signal),20);
        set(h,'LineColor','none')
        cmap=[0.182352945208550,0.182352945208550,0.790196120738983;0.174424752593040,0.211379334330559,0.802652835845947;...
            0.166161164641380,0.242508038878441,0.815109610557556;0.157562181353569,0.275798201560974,0.827566325664520;...
            0.148627817630768,0.311309069395065,0.840023100376129;0.139358073472977,0.349099755287170,0.852479815483093;...
            0.129752933979034,0.389229506254196,0.864936590194702;0.119812399148941,0.431757479906082,0.877393305301666;...
            0.109536483883858,0.476742863655090,0.889850080013275;0.0989251807332039,0.524244844913483,0.902306795120239;...
            0.0879784896969795,0.574322640895844,0.914763569831848;0.0766964107751846,0.627035379409790,0.927220284938812;...
            0.0650789439678192,0.682442307472229,0.939677059650421;0.0531260855495930,0.740602552890778,0.952133774757385;...
            0.0408378429710865,0.801575362682343,0.964590549468994;0.0282142106443644,0.865419864654541,0.977047264575958;...
            0.0152551913633943,0.932195305824280,0.989504039287567;0.00196078442968428,1,1;0.00196078442968428,1,0.911051690578461;...
            0.00196078442968428,1,0.820142626762390;0.00196078442968428,1,0.729233503341675;0.00196078442968428,1,0.638324439525604;...
            0.00196078442968428,1,0.547415316104889;0.00196078442968428,1,0.456506252288818;0.00196078442968428,1,0.365597158670425;...
            0.00196078442968428,1,0.274688065052032;0.00196078442968428,1,0.183778971433640;0.00196078442968428,1,0.0928698778152466;...
            0.00196078442968428,1,0.00196078442968428;0.0928698778152466,1,0.00196078442968428;0.183778971433640,1,0.00196078442968428;...
            0.274688065052032,1,0.00196078442968428;0.365597158670425,1,0.00196078442968428;0.456506252288818,1,0.00196078442968428;...
            0.547415316104889,1,0.00196078442968428;0.638324439525604,1,0.00196078442968428;0.729233503341675,1,0.00196078442968428;...
            0.820142626762390,1,0.00196078442968428;0.911051690578461,1,0.00196078442968428;1,1,0.00196078442968428;1,0.939460813999176,0.00196078442968428;...
            1,0.876960813999176,0.00196078442968428;1,0.814460813999176,0.00196078442968428;1,0.751960813999176,0.00196078442968428;...
            1,0.689460813999176,0.00196078442968428;1,0.626960813999176,0.00196078442968428;1,0.564460813999176,0.00196078442968428;...
            1,0.501960813999176,0.00196078442968428;1,0.439460784196854,0.00196078442968428;1,0.376960784196854,0.00196078442968428;...
            1,0.314460784196854,0.00196078442968428;1,0.251960784196854,0.00196078442968428;1,0.189460784196854,0.00196078442968428;...
            1,0.126960784196854,0.00196078442968428;1,0.0644607841968536,0.00196078442968428;1,0.00196078442968428,0.00196078442968428;...
            0.939460813999176,0.00196078442968428,0.00196078442968428;0.876960813999176,0.00196078442968428,0.00196078442968428;...
            0.814460813999176,0.00196078442968428,0.00196078442968428;0.751960813999176,0.00196078442968428,0.00196078442968428;...
            0.689460813999176,0.00196078442968428,0.00196078442968428;0.626960813999176,0.00196078442968428,0.00196078442968428;...
            0.564460813999176,0.00196078442968428,0.00196078442968428;0.501960813999176,0.00196078442968428,0.00196078442968428];
        colormap(cmap);
        xlabel('Emission (nm)');  ylabel('Excitation (nm)'); 
        set(gca,'fontweight','bold','fontsize',15)
        set(gca,'PlotBoxAspectRatio',[1 1 1])
        set(gca,'XTick',[200,300,400,500,600,700,800]);
        set(gca,'YTick',[200,250,300,350,400,450,500,550,600,650,700,750,800]);
        if Peak_Labels 
            PeakLabels_color='k';
            PeakLabels_fontsize_simple=20;
            text(305,275,'B','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)
            text(340,275,'T','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)    
            text(370,280,'N','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)    
            text(395,300,'M','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)    
            text(450,255,'A','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)
            text(445,342.5,'C','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)
            text(509,390,'D','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)
            %text(521,455,'E','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)    
            %text(660,398,'P','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)

            % Can disable all lines above and enable all bottom ones for more detailed peak labels
%             PeakLabels_fontsize_advanced=15;
%             text(305,230,'B_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)        
%             text(305,275,'B_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(340,230,'T_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(340,275,'T_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             text(370,280,'N','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             text(375,240,'M_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(395,300,'M_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             text(430,260,'C_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(445,342.5,'C_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             text(487,250,'C^{+}_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(487,402.5,'C^{+}_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)  
%             text(509,390,'D','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             %text(521,455,'E','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             %text(660,398,'P','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             %text(312.5,230,'H','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
        end
            handle_colorbar=colorbar;
            set(get(handle_colorbar,'title'),'string','RU','fontsize',13);
end

% Plotting EEM with x-axis = excitation
function plotEEM_SwitchedAxis(tit,Signal,Em,Ex)
global  Peak_Labels
        title(tit,'Interpreter','none','fontsize',13,'Position',[median(Ex),max(Em)+10,0]);
        [~,h]=contourf(Ex,Em,squeeze(Signal'),20);
        set(h,'LineColor','none')
        cmap=[0.182352945208550,0.182352945208550,0.790196120738983;0.174424752593040,0.211379334330559,0.802652835845947;...
            0.166161164641380,0.242508038878441,0.815109610557556;0.157562181353569,0.275798201560974,0.827566325664520;...
            0.148627817630768,0.311309069395065,0.840023100376129;0.139358073472977,0.349099755287170,0.852479815483093;...
            0.129752933979034,0.389229506254196,0.864936590194702;0.119812399148941,0.431757479906082,0.877393305301666;...
            0.109536483883858,0.476742863655090,0.889850080013275;0.0989251807332039,0.524244844913483,0.902306795120239;...
            0.0879784896969795,0.574322640895844,0.914763569831848;0.0766964107751846,0.627035379409790,0.927220284938812;...
            0.0650789439678192,0.682442307472229,0.939677059650421;0.0531260855495930,0.740602552890778,0.952133774757385;...
            0.0408378429710865,0.801575362682343,0.964590549468994;0.0282142106443644,0.865419864654541,0.977047264575958;...
            0.0152551913633943,0.932195305824280,0.989504039287567;0.00196078442968428,1,1;0.00196078442968428,1,0.911051690578461;...
            0.00196078442968428,1,0.820142626762390;0.00196078442968428,1,0.729233503341675;0.00196078442968428,1,0.638324439525604;...
            0.00196078442968428,1,0.547415316104889;0.00196078442968428,1,0.456506252288818;0.00196078442968428,1,0.365597158670425;...
            0.00196078442968428,1,0.274688065052032;0.00196078442968428,1,0.183778971433640;0.00196078442968428,1,0.0928698778152466;...
            0.00196078442968428,1,0.00196078442968428;0.0928698778152466,1,0.00196078442968428;0.183778971433640,1,0.00196078442968428;...
            0.274688065052032,1,0.00196078442968428;0.365597158670425,1,0.00196078442968428;0.456506252288818,1,0.00196078442968428;...
            0.547415316104889,1,0.00196078442968428;0.638324439525604,1,0.00196078442968428;0.729233503341675,1,0.00196078442968428;...
            0.820142626762390,1,0.00196078442968428;0.911051690578461,1,0.00196078442968428;1,1,0.00196078442968428;1,0.939460813999176,0.00196078442968428;...
            1,0.876960813999176,0.00196078442968428;1,0.814460813999176,0.00196078442968428;1,0.751960813999176,0.00196078442968428;...
            1,0.689460813999176,0.00196078442968428;1,0.626960813999176,0.00196078442968428;1,0.564460813999176,0.00196078442968428;...
            1,0.501960813999176,0.00196078442968428;1,0.439460784196854,0.00196078442968428;1,0.376960784196854,0.00196078442968428;...
            1,0.314460784196854,0.00196078442968428;1,0.251960784196854,0.00196078442968428;1,0.189460784196854,0.00196078442968428;...
            1,0.126960784196854,0.00196078442968428;1,0.0644607841968536,0.00196078442968428;1,0.00196078442968428,0.00196078442968428;...
            0.939460813999176,0.00196078442968428,0.00196078442968428;0.876960813999176,0.00196078442968428,0.00196078442968428;...
            0.814460813999176,0.00196078442968428,0.00196078442968428;0.751960813999176,0.00196078442968428,0.00196078442968428;...
            0.689460813999176,0.00196078442968428,0.00196078442968428;0.626960813999176,0.00196078442968428,0.00196078442968428;...
            0.564460813999176,0.00196078442968428,0.00196078442968428;0.501960813999176,0.00196078442968428,0.00196078442968428];
        colormap(cmap);
        xlabel('Excitation (nm)');  ylabel('Emission (nm)'); 
        set(gca,'fontweight','bold','fontsize',15)
        set(gca,'PlotBoxAspectRatio',[1 1 1])
        set(gca,'YTick',[200,300,400,500,600,700,800]);
        set(gca,'XTick',[200,250,300,350,400,450,500,550,600,650,700,750,800]);
        if Peak_Labels 
            PeakLabels_color='k';
            PeakLabels_fontsize_simple=20;
            text(275,305,'B','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)
            text(275,340,'T','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)    
            text(280,370,'N','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)    
            text(300,395,'M','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)    
            text(255,450,'A','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)
            text(342.5,445,'C','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)
            text(390,509,'D','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)
            %text(455,521,'E','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)    
            %text(398,660,'P','fontweight','bold','fontsize',PeakLabels_fontsize_simple,'Color',PeakLabels_color)

            % Can disable all lines above and enable all bottom ones for more detailed peak labels
%             PeakLabels_fontsize_advanced=15;
%             text(230,305,'B_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(275,305,'B_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(230,340,'T_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(275,340,'T_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             text(280,370,'N','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             text(240,375,'M_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(300,395,'M_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             text(260,430,'C_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(342.5,445,'C_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             text(250,487,'C^{+}_1','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(402.5,487,'C^{+}_2','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(390,509,'D','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(455,521,'E','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             text(660,398,'P','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             text(230,312.5,'H','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
        end
            handle_colorbar=colorbar;
            set(get(handle_colorbar,'title'),'string','RU','fontsize',13);
end