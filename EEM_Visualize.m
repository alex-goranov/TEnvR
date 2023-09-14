function EEM_Visualize(filename,varargin)
%% Description

% This code takes one, two or three processed EEM spectra and plots them. 
% The spectra must be located in .csv files (comma-separated values). 

% Examples:
% EEM_Visualize('Sample 1_ref_Final.csv')
% EEM_Visualize('Sample 1_ref_Final.csv','Sample 2_ref_Final.csv')
% EEM_Visualize('Sample 1_ref_Final.csv','Sample 2_ref_Final.csv','Sample 3_ref_Final.csv')

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

%% Configuration
clearvars -except filename varargin
close all

global Same_Intensity_Scale Peak_Labels bottom top % Do not edit this line and do not delete bottom/top!

% Peak Labels
Peak_Labels=true; % choose true or false

% EEM orientation: choose if you want the excitation to be as x-axis
Excitation_Xaxis=false; % choose true if you want excitation to be as x-exis; false = y-axis

% Same_Intensity_Scale: if you want all displayed EEMs to be with the same intensity scales
Same_Intensity_Scale=false; % choose true or false

% Normalization: choose  if you want to normalize the EEMs to their total spectral intensity
Normalization=false; % choose true or false

%% Option 1: Plotting one spectrum

if nargin==1 % EEM_Visualize(filename)
    Data=csvread(filename); % Import data
    Ex=Data(1,:); Em=Data(:,1); % Extract excitation and emission ranges
    Ex=Ex(2:end); Em=Em(2:end); 
    Int=Data; Int(1,:)=[]; Int(:,1)=[]; % Extract intensity data
    Int=transpose(Int);
    
    % Make sure the dimensions are consistent
    Size_Em=size(Em);
    Size_Ex=size(Ex);
    Size_Int=size(Int);
    if Size_Em(1)~= Size_Int(2)
        error('Error! Check the dimensions of the emission of sample 1')
    end
    if Size_Ex(2)~= Size_Int(1)
        error('Error! Check the dimensions of the excitation of sample 1')
    end
    
    if Normalization
        Int=Int/nansum(nansum(Int));
    end
    
    % In case Same_Intensity_Scale was mistakenly turned on:
    Same_Intensity_Scale=false;
    
    figure(1);
    figurename=char(filename(1:end-4)); % Figure name is sample as sample name without ".csv"
        hold on
        if Excitation_Xaxis 
            plotEEM_SwitchedAxis(char(filename(1:end-4)),Int,Em,Ex)
        else
            plotEEM(char(filename(1:end-4)),Int,Em,Ex) 
        end
        hold off
        set(gcf,'Units','inches')
        set(gcf,'InnerPosition',[1,1,7,7]);
        set(gcf,'OuterPosition',[1,1,8,8]);
        print(gcf,['EEM_' figurename '.tiff'],'-dtiff','-r500');
        disp(['Finished plotting ' char(filename) '!'])

%% Option 2: Plotting two spectra 

elseif nargin==2 % EEM_Visualize(filename,filename2)
    filename1=filename;
    filename2=varargin{1};
    Data1=csvread(filename1); % Import data
    Data2=csvread(filename2); % Import data
    Ex1=Data1(1,:); Em1=Data1(:,1); Ex2=Data2(1,:); Em2=Data2(:,1);  % Extract excitation and emission ranges
    Ex1=Ex1(2:end); Em1=Em1(2:end); Ex2=Ex2(2:end); Em2=Em2(2:end);  % Extract excitation and emission ranges
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

    if Normalization
        Int1=Int1/nansum(nansum(Int1));
        Int2=Int2/nansum(nansum(Int2));
    end
    
    if Same_Intensity_Scale
        top=max(max([Int1,Int2]));
        bottom=min(min([Int1,Int2]));
        if bottom < 0
            error('Error! The EEM contains a negative number!')
        end
        caxis manual
        caxis([bottom top]);
    end
    
    % Plotting
    figure(1);
    figurename=[char(filename1(1:end-4)) '_' char(filename2(1:end-4))]; % Figure name is contain both sample names without ".csv"
    f1=subplot(1,2,1);
        hold on
        if Excitation_Xaxis 
            plotEEM_SwitchedAxis(char(filename1(1:end-4)),Int1,Em1,Ex1)
        else
            plotEEM(char(filename1(1:end-4)),Int1,Em1,Ex1)
        end
        hold off
    f2=subplot(1,2,2);
        hold on
        if Excitation_Xaxis 
            plotEEM_SwitchedAxis(char(filename2(1:end-4)),Int2,Em2,Ex2)
        else
            plotEEM(char(filename2(1:end-4)),Int2,Em2,Ex2)
        end
        hold off
        set(gcf,'Units','inches')
        set(gcf,'OuterPosition',[-0.104166666666667,3.864583333333333,14.864583333333332,7.385416666666667]);
        set(f1,'Position',[0.06,0.13,0.34,0.78]);
        set(f2,'Position',[0.57,0.13,0.34,0.78]);
        print(gcf,['EEM_' figurename '.tiff'],'-dtiff','-r500');
        disp(['Finished plotting ' char(filename1) ' and ' char(filename2) '!'])

%% Option 3: Plotting 3 spectra

elseif nargin==3 % EEM_Visualize(filename,filename2,filename3)
    filename1=filename;
    filename2=varargin{1};
    filename3=varargin{2};
    Data1=csvread(filename1); % Import data
    Data2=csvread(filename2); % Import data
    Data3=csvread(filename3); % Import data
    Ex1=Data1(1,:); Em1=Data1(:,1); Ex2=Data2(1,:); Em2=Data2(:,1); Ex3=Data3(1,:); Em3=Data3(:,1); % Extract excitation and emission ranges
    Ex1=Ex1(2:end); Em1=Em1(2:end); Ex2=Ex2(2:end); Em2=Em2(2:end); Ex3=Ex3(2:end); Em3=Em3(2:end); % Extract excitation and emission ranges
    Int1=Data1; Int2=Data2; Int3=Data3; Int1(1,:)=[]; Int1(:,1)=[]; Int2(1,:)=[]; Int2(:,1)=[]; Int3(1,:)=[]; Int3(:,1)=[]; % Extract intensity data
    Int1=transpose(Int1); Int2=transpose(Int2); Int3=transpose(Int3);

    % Make sure the dimensions are consistent
    Size_Em1=size(Em1);
    Size_Em2=size(Em2);
    Size_Em3=size(Em3);
    Size_Ex1=size(Ex1);
    Size_Ex2=size(Ex2);
    Size_Ex3=size(Ex3);
    Size_Int1=size(Int1);
    Size_Int2=size(Int2);
    Size_Int3=size(Int3);

    if Size_Em1(1)~= Size_Int1(2)
        error('Error! Check the dimensions of the emission of sample 1')
    end
    if Size_Em2(1)~= Size_Int2(2)
        error('Error! Check the dimensions of the emission of sample 2')
    end
    if Size_Em3(1)~= Size_Int3(2)
        error('Error! Check the dimensions of the emission of sample 3')
    end
    if Size_Ex1(2)~= Size_Int1(1)
        error('Error! Check the dimensions of the excitation of sample 1')
    end
    if Size_Ex2(2)~= Size_Int2(1)
        error('Error! Check the dimensions of the excitation of sample 2')
    end
    if Size_Ex3(2)~= Size_Int3(1)
        error('Error! Check the dimensions of the excitation of sample 3')
    end

    if Normalization
        Int1=Int1/nansum(nansum(Int1));
        Int2=Int2/nansum(nansum(Int2));
        Int3=Int3/nansum(nansum(Int3));
    end
    
    if Same_Intensity_Scale
            top=max(max([Int1,Int2,Int3]));
            bottom=min(min([Int1,Int2,Int3]));
            if bottom < 0
                error('Error! The EEM contains a negative number!')
            end
        caxis manual
        caxis([bottom top]);
    end
       
    % Plotting
    figure(1);
    figurename=[char(filename1(1:end-4)) '_' char(filename2(1:end-4)) '_' char(filename3(1:end-4))]; % Figure name is contain all sample names without ".csv"
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
            hold on
            if Excitation_Xaxis
                plotEEM_SwitchedAxis(char(filename3(1:end-4)),Int3,Em3,Ex3)
            else
                plotEEM(char(filename3(1:end-4)),Int3,Em3,Ex3)
            end
        hold off
        screen=get(0,'screensize');
        set(gcf,'OuterPosition',[0,screen(4)/2,screen(3),screen(4)/2]);
        print(gcf,['EEM_' figurename],'-dtiff','-r500');
        disp(['Finished plotting ' char(filename1) ', ' char(filename2) ', and ' char(filename3) '!'])
else
    error('Error! You have entered too many variables for this function!')
end
end

%% Internal Functions

% Plotting EEM with x-axis = emission
function plotEEM(tit,Signal,Em,Ex)
global  Same_Intensity_Scale Peak_Labels bottom top 
        if Same_Intensity_Scale
            caxis manual
            caxis([bottom top]);
        end
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
        else
            ;
        end
        handle_colorbar=colorbar;
        if ~contains(tit,'comp')
            set(get(handle_colorbar,'title'),'string','RU','fontsize',13);
        end
end

% Plotting EEM with x-axis = excitation
function plotEEM_SwitchedAxis(tit,Signal,Em,Ex)
global  Same_Intensity_Scale Peak_Labels bottom top 
        if Same_Intensity_Scale
            caxis manual
            caxis([bottom top]);
        end      
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
%             %text(455,521,'E','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)    
%             %text(660,398,'P','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
%             %text(230,312.5,'H','fontweight','bold','fontsize',PeakLabels_fontsize_advanced,'Color',PeakLabels_color)
        else
            ;
        end
        handle_colorbar=colorbar;
        if ~contains(tit,'comp')
            set(get(handle_colorbar,'title'),'string','RU','fontsize',13);
        end
end