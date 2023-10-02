%% Description

% This code loads multiple formula lists, identifies and identifies common
% and unique formulas among each sample. Creates an alignment matrix. 
% Prepares data for multivariate statistics. 

% The alignment algorithm was modified by Krista Longnecker and Elizabeth Kujawinski 
% (Woods Hole Oceanographic Institution) from code originally published in 
% Mantini et al. (2007): 

% Mantini,Petrucci,Pieragostino, Del Boccio, Di Nicola, Di Ilio, Federici, Sacchetta, Comani and Urbani (2007). 
% "LIMPIC: a computational method for the separation of protein MALDI-TOF-MS signals from noise." BMC Bioinformatics 8: 101.

% Example: 
% Run the code by typing FTMS_AlignmentFormulas in the command windwow.

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

%% Import FTMS data
clear
close all
tic
[data, filename] = FTMSListsLoad; 

for i=1:size(data,2)
    M(i)=size(data{i},2);
end

if std(M) ~=0
    error('Error! There is a file of inconsistent format in the Current Folder! Remove it.')
end

global config

config=FTMS_ConfigurationToolbox;

%% Align peaks and Export to Excel
[Peaks, Intensity, AllFormulas] = SpectralAlignment(data,0.1,config.Column_ExactMass,config.Column_Magnitude);

matrix_alignedpeaks=[Peaks Intensity];

exportname='FTMS Alignment Matrix.xlsx';
warning('off','MATLAB:xlswrite:AddSheet')
xlswrite(exportname,['ExactMass' filename],1,'A1');
xlswrite(exportname,matrix_alignedpeaks,1,'A2'); 
 
%% Align Molecular parameters with peaks, Trim, and Export to Excel

Search_Masses  = AllFormulas(:,config.Column_ExactMass);
Search_M  = zeros(1, size(Peaks,1));

   for Search_iA = 1:size(Peaks,1)
      Search_dist = abs(Peaks(Search_iA) - Search_Masses);
      SearchInd = find(ismember(Search_dist,min(Search_dist),'rows'));
      
      if size(SearchInd,1) ~=1
          SearchInd=SearchInd(1,1);
      end

      if ~isempty(SearchInd)
         Search_M(Search_iA) = SearchInd;
      end
   end
   
   AlignedFormulas=[];
   for Search_i=1:size(Search_M,2)
       AlignedFormulas(Search_i,:)=AllFormulas(Search_M(Search_i),:);
   end
   
% Trim useless parameters 
AlignedFormulas(:,[config.Column_mz,config.Column_Magnitude,config.Column_SN,config.Column_EstimC,config.Column_Type,config.Column_ExactMass,config.Column_error])= [];
   
matrix_alignedpeakswithformulas=[Peaks, Intensity, AlignedFormulas];

xlswrite(exportname,['ExactMass' filename 'C' 'H' 'O' 'N' 'S' 'P' 'E' 'O/C' 'H/C' 'DBE' 'DBE/C' 'AImod'],'Formulas','A1');
xlswrite(exportname,matrix_alignedpeakswithformulas,'Formulas','A2');

% Trimming the matrix based on user-defined cut-offs (e.g., m/z 300-800)
indeces_mz = matrix_alignedpeakswithformulas(:,1)>=config.MZ_low & matrix_alignedpeakswithformulas(:,1)<=config.MZ_high;
matrix_alignedALL_trimmed=matrix_alignedpeakswithformulas(indeces_mz,:);

% Trimming the matrix based formulas found in minimum number of samples
FormulasFoundInSamples=sum(matrix_alignedALL_trimmed(:,2:end-12)>0,2);
indeces_NumSamples=FormulasFoundInSamples >= config.MinSamples;
matrix_alignedALL2_trimmed=matrix_alignedALL_trimmed(indeces_NumSamples,:);

xlswrite(exportname,['ExactMass' filename 'C' 'H' 'O' 'N' 'S' 'P' 'E' 'O/C' 'H/C' 'DBE' 'DBE/C' 'AImod'],'Trimmed','A1');
xlswrite(exportname,matrix_alignedALL2_trimmed,'Trimmed','A2');

%% Normalize & Export to Excel

Intensity=matrix_alignedALL2_trimmed(:,2:end-12);

for i=1:size(Intensity,2)
    IntValues=Intensity(:,i);
    Totalint=sum(IntValues);
    IntValues=IntValues./Totalint;
    Intensity_norm(:,i)=IntValues;
end

if config.PresenceAbsence
Intensity_norm(Intensity_norm > 0) = 1;
Intensity_norm=Intensity_norm./sum(Intensity_norm,1);
end

matrix_normalized=[matrix_alignedALL2_trimmed(:,1), Intensity_norm, matrix_alignedALL2_trimmed(:,end-11:end)];
xlswrite(exportname,['ExactMass' filename 'C' 'H' 'O' 'N' 'S' 'P' 'E' 'O/C' 'H/C' 'DBE' 'DBE/C' 'AImod'],'Normalized','A1');
xlswrite(exportname,matrix_normalized,'Normalized','A2');

%% Determine Common Formulas & Export to Excel

matrix_common=matrix_normalized(:,1:end-12);
matrix_common(matrix_common == 0) = NaN;
matrix_common=[matrix_common, matrix_alignedALL2_trimmed(:,end-11:end)];
matrix_common(any(isnan(matrix_common), 2), :) = [];
if ~isempty(matrix_common)
    xlswrite(exportname,['ExactMass' filename 'C' 'H' 'O' 'N' 'S' 'P' 'E' 'O/C' 'H/C' 'DBE' 'DBE/C' 'AImod'],'Normalized Common','A1'); 
    xlswrite(exportname,matrix_common,'Normalized Common','A2'); 
end
disp(['Finished aligning formulas (' num2str(toc) ' seconds)'])

%% Internal Functions

% Function for loading mass lists
function [data, filename] = FTMSListsLoad
D=dir('*_Final.xlsx'); %Adjust for different file extensions
D=natsortfiles(D);

if isempty(D)
    error('Error! Something is wrong because no files were found!')
end

for a=1:size(D,1)
    fname = [D(a).name];
    
    data{a}=xlsread(fname);
    data{a}=sortrows(data{a},1);
    filename{a}=[D(a).name(1:end-5)];
end

end

% Alignment algorhitm
function [Peaks, Intensity, data_allformulas]=SpectralAlignment(dataIn,ClusteringTolerance,ColumnMass,ColumnInt)

global config

data =[];
data_allstacked=[];
for a=1:size(dataIn,2)
    peak_mass = dataIn{a}(:,ColumnMass);
    peak_intensity = dataIn{a}(:,ColumnInt);
    ntp=size(peak_intensity,1);
    data=[data ; [peak_mass' ; peak_intensity'; a*ones(1,ntp)]']; % compiles all mass lists together
    data_allstacked = [data_allstacked; dataIn{a}];
    clear ntp peak_mass peak_intensity 
end

data_allformulas=data_allstacked;
[~,ia,~] = unique(data_allformulas(:,config.Column_ExactMass),'rows');
data_allformulas=data_allformulas(ia,:);
data_allformulas=sortrows(data_allformulas,config.Column_ExactMass);
 
nspectra = size(unique(data(:,3)),1);
% creation of data matrices for peak clustering
[dataord(:,1),order]=sort(data(:,1));
dataord(:,2:3)=data(order,2:3);

Peaks=dataord(:,1)';
npoints=size(Peaks,2);
lista_nel=ones(1,npoints);
Intensity=zeros(nspectra,npoints);
for i=1:npoints
    pat=dataord(i,3);
    amp=dataord(i,2);
    Intensity(pat,i)=amp;
end

% iterative peak clustering
i=2;

while i < size(Peaks,2)
    
    value_b=Peaks(i-1);
    value=Peaks(i);
    value_p=Peaks(i+1);
    db=value-value_b;
    dp=value_p-value;
    thres=value*ClusteringTolerance/1000000;
    pat=find(Intensity(:,i) > 0); 
    logicval=true;
    
    % looking at the peak on the lower end
    for k=1:size(pat,1)
        logic =(Intensity(pat(k), i-1 ) == 0) ;
        logicval=logicval & logic;
    end

    if logicval == false
        db = inf;
    end
    logicval=true;
    
    % Now do the same for the peak on the upper end
    for k=1:size(pat,1)
        logic =(Intensity(pat(k), i+1 ) == 0) ;
        logicval=logicval & logic;
    end
    
    if logicval == false
        dp = inf;
    end
    
    if db <= dp && db < thres    
      % peak1 and peak2 are closer together:
      % Take the average m/z of peak1 and peak2 and put that into Peaks
      Peaks(i-1)=(lista_nel(i)*value+lista_nel(i-1)*value_b)/(lista_nel(i-1)+lista_nel(i));
      lista_nel(i-1)=lista_nel(i-1)+lista_nel(i); 
      
      % Adding the intensity of peak1 and peak2
      Intensity(:,i-1)=Intensity(:,i-1)+Intensity(:,i);
      
      % and delete the index
      Peaks(i)=[];
      lista_nel(i)=[];
      Intensity(:,i)=[];
    elseif db > dp && dp < thres 
      % peak2 and peak3 are closer together (bc distance between peak1 and peak2 is greater):
      % Take the average m/z of peak2 and peak3 and put that into Peaks
      Peaks(i+1)=(lista_nel(i)*value+lista_nel(i+1)*value_p)/(lista_nel(i+1)+lista_nel(i));
      lista_nel(i+1)=lista_nel(i+1)+lista_nel(i); 
      
      % Adding the intensity of peak2 and peak3
      Intensity(:,i+1)=Intensity(:,i+1)+Intensity(:,i);
      
      % and delete the index
      Peaks(i)=[];
      lista_nel(i)=[];
      Intensity(:,i)=[];  
    else
        % the distance bw the two peaks is greater than the threshold - move on to the next set
        i=i+1;
    end
end

lista_min=Peaks*(1-ClusteringTolerance/1000000); 
lista_max=Peaks*(1+ClusteringTolerance/1000000);

for i=2:size(Peaks,2)
    if lista_min(i) < lista_max(i-1)
        lista_min(i)=mean(Peaks(i-1:i));
        lista_max(i-1)=lista_min(i);
    end
end

Peaks = Peaks';
Intensity = Intensity';
end