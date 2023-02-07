function FTMS_RefinementFormulas(filename,config)
%% Description

% This code takes an list of assigned formulas and:
    % Stage 1: Import Data & reformat  
    % Stage 2: H-correction, Reformatting, and initial filtering using Elemental Constrains
    % Stage 3: Formula refinement using Isotopic filter and selection of Unique formulas
    % Stage 4: Formula refinement using KMD Filter
    % Stage 5: Formula refinement using Compositional filter
    % Stage 6: Formula refinement using Error filter
    % Stage 7: Final refinement (removal of residual ambiguous assignments)
    % Stage 8: Quality Control (produces a figure)

% Examples:
% FTMS_RefinementFormulas('Sample1_Refinement_F.txt',FTMS_ConfigurationFile)

%% Copyright and License Notice: 

% Copyright © 2022 Old Dominion University Research Foundation, Norfolk VA, USA
% All rights reserved.

% This file is part of the Toolbox for Environmental Research (TEnvR). Please cite the toolbox as follows: 
% Goranov, A. I., Sleighter, R. L., Yordanov, D. A., and Hatcher, P. (2023): 
% TEnvR: MATLAB-Based Toolbox for Environmental Research, Journal TBD, doi: XXXXXXXXXXX.

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

%% Stage 1: Import Data & Reorganize file from Molecular Formula Calculators 
tic

fid=fopen(filename);    % open data file
str=fgetl(fid);         % get first line in data file

if str(1)=='O'          % Loads and reformats the file from FTMS_FormulaAssignment
    Data_Stage1=readtable(filename,'VariableNamingRule','preserve');
    fclose(fid);
    Data_Stage1=table2array(Data_Stage1);
    for i=1:size(Data_Stage1,1)
        if i == 1
            Index(i)=1;
            Previous_Index=Index(i);
        elseif i < size(Data_Stage1,1) && i~=1 
            if Data_Stage1(i,2) ~= Data_Stage1(i-1,2) 
                Index(i)=Previous_Index+1;
                Previous_Index=Index(i);
            elseif Data_Stage1(i,2) == Data_Stage1(i-1,2) 
                Index(i)=Previous_Index;
            end
        elseif i == size(Data_Stage1,1)
            if Data_Stage1(i,2) ~= Data_Stage1(i-1,2)
                Index(i)=Previous_Index+1;
            elseif Data_Stage1(i,2) == Data_Stage1(i-1,2)
                Index(i)=Previous_Index;
            end
        end
    end
Index=Index';    
Data_Stage1=[Index,Data_Stage1(:,2),zeros(size(Data_Stage1,1),1)+config.MassAccuracy,Data_Stage1(:,3:14)];

else % Loads and reformats the file from the Molecular Formula Calculator (obsolete)

n=0;
i=1;
while ischar(str)               % end of line str outputs -1 and is not char
    if ~isempty(str)            % empty strings are line with no data
        if str(1)=='P'          % header lines begin with 'P' in Peak report
            hdr=retheader(str); % use the internal retheader function to return numerical values from header line
            n=n+1;
        elseif str(1)=='E'
            error('ERROR! Number of iterations exceeded! Re-assign formulas with MFC and add more zeros to "Maximum Iterations"!')
        else         
            dat=retdata(str);    % use redata subfuction to return numerical values from data lines
            data=[n hdr' dat];   % concatenate index, header info, and data
            Data_Stage1(i,:)=data;
            i=i+1;
        end
    end  
    str=fgetl(fid);             % get next line in data file
end
fclose(fid);

Data_Stage1(isnan(Data_Stage1))=0; % substitute all NaNs with zeros and 
end

% At this point the data from either input (MFC vs MATLAB) is formatted in a consistent way

Data_Stage1=sortrows(Data_Stage1,2); % Sort based on m/z, column 2

% Import S/N and C-number data from the Refinement file
RefinementData= xlsread([filename(1:end-6) '.xlsx'], 'DataRefined');        % Import
RefinementQC= xlsread([filename(1:end-6) '.xlsx'], 'Sheet1','AH3:AJ13');    % Import
RefinementData=sortrows(RefinementData,1); % Sort based on m/z, column 1

% Align S/N and C-number data from the Refinement file into the formulas file
i = 1; j = 1;
while i <= size(RefinementData,1) && j <= size(Data_Stage1, 1)
    mz_BeforeAssignment = RefinementData(i, 1);
    mz_AfterAssignment = Data_Stage1(j, 2);

    if intersect(mz_BeforeAssignment,mz_AfterAssignment) > 0
        Data_Stage1_Expanded(j,:)=[Data_Stage1(j,:),RefinementData(i,3:4)];
        j = j+1;
    elseif mz_BeforeAssignment < mz_AfterAssignment
        i = i+1;
    elseif mz_BeforeAssignment > mz_AfterAssignment
        j = j+1;
    end
end

% Export
warning('off','MATLAB:xlswrite:AddSheet');
xlswrite([filename(1:end-6) '_Processing.xlsx'],{'Index' 'm/z' 'Mass Accuracy (ppm)' 'Magnitude' 'C' 'Hion' 'N' 'O' 'S' 'P' 'E' 'K' 'Na' 'm/z corr' 'Assignment Error (ppm)','S/N','Estim C#'},'Sheet1','A1');
xlswrite([filename(1:end-6) '_Processing.xlsx'],Data_Stage1_Expanded,'Sheet1','A2');

xlswrite([filename(1:end-6) '_Processing.xlsx'],{'A-O is the reformatted output from the Formula Assignment (either MFC or MATLAB)'},'Sheet1','S2');
xlswrite([filename(1:end-6) '_Processing.xlsx'],{'P-Q is data from the _Refinement.xlsx file!'},'Sheet1','S3');
xlswrite([filename(1:end-6) '_Processing.xlsx'],{'Hion is number of H atoms in the singly charged ion!'},'Sheet1','S4');
xlswrite([filename(1:end-6) '_Processing.xlsx'],{'This is the initial formatting of the data, data hereafter will be formatted differently!'},'Sheet1','S5');

%% Stage 2: H-correction, Reformatting, and Initial Filtering using Elemental Constraints

clearvars -except filename config RefinementData RefinementQC Data_Stage1_Expanded 

Data_Stage1_Expanded=sortrows(Data_Stage1_Expanded,2);

% Pick the imported data apart into individual arrays
Index=Data_Stage1_Expanded(:,1);
mz=Data_Stage1_Expanded(:,2);
%AssignmentAccuracy=Data_Stage1_Expanded(:,3); % This is not used in further calculations.
Magnitude=Data_Stage1_Expanded(:,4);
C=Data_Stage1_Expanded(:,5);
Hion=Data_Stage1_Expanded(:,6);
N=Data_Stage1_Expanded(:,7);
O=Data_Stage1_Expanded(:,8);
S=Data_Stage1_Expanded(:,9);
P=Data_Stage1_Expanded(:,10);
E=Data_Stage1_Expanded(:,11);
K=Data_Stage1_Expanded(:,12);
Na=Data_Stage1_Expanded(:,13);
%ExactMass_ion=Data_Stage1_Expanded(:,14); % This is not used in further calculations.
e=Data_Stage1_Expanded(:,15);
e=round(e,1);
SN=Data_Stage1_Expanded(:,16);
EstimC=Data_Stage1_Expanded(:,17);

% Correct the number of H based on the ionization mode (i.e., convert ion formulas into molecular formulas)
if strcmpi(config.Ionization_Mode,'negative')
    H=Hion+1;
elseif strcmpi(config.Ionization_Mode,'positive')
    H=Hion-1+Na+K; % Keep in mind that there are some peaks which were assigned with both Na and K!!!
else
    error('Error! Need to select the correct ionization mode (positive or negative) in the configuration file!')
end

% Calculate the ExactMass of the molecule using the mass of each element
ExactMass=C*config.Mass_12C+H*config.Mass_1H+N*config.Mass_14N+O*config.Mass_16O+S*config.Mass_32S+P*config.Mass_31P+E*config.Mass_Heteroelement;
ExactMass=round(ExactMass,config.Precision+1); % adds an uncertainty of 1 decimal to the defined precision

% Calculate other important parameters (O/C, H/C, DBE, DBE/C, AImod)
OC=O./C;
HC=H./C;
    
if config.Heteroelement_Halogen && ~config.Heteroelement_N15        % E = 35Cl
    DBE=1+(((2*C)-H+N+P-E)/2);
    AImod=(1+C-(0.5.*O)-S-(0.5.*(N+P+H+E)))./(C-(0.5.*O)-N-S-P);
elseif ~config.Heteroelement_Halogen && config.Heteroelement_N15    % E = 15N
    DBE=1+(((2*C)-H+N+P+E)/2);
    AImod=(1+C-(0.5.*O)-S-(0.5.*(N+P+H+E)))./(C-(0.5.*O)-N-S-P-E);
else %Calculate using only CHONSP elements
    DBE=1+(((2*C)-H+N+P)/2);                                            
    AImod=(1+C-(0.5.*O)-S-(0.5.*(N+P+H)))./(C-(0.5.*O)-N-S-P);
end

AImod(AImod <= 0)=0; 
DBEC=DBE./C;

% Filtering using Elemental Constraints

NumberCases      = 18; % Number of cases must match number of rejection criteria below!
Rejected  = zeros([size(Index,1),NumberCases]);

Rejected(:,1) = OC < config.O_to_C_min;
Rejected(:,2) = OC > config.O_to_C_max;
Rejected(:,3) = O-C > config.O_to_C_delta_max;

Rejected(:,4) = HC < config.H_to_C_min;
Rejected(:,5) = HC > config.H_to_C_max;

Rejected(:,6) = N./C < config.N_to_C_min;
Rejected(:,7) = N./C > config.N_to_C_max;

Rejected(:,8) = S./C < config.S_to_C_min;
Rejected(:,9) = S./C > config.S_to_C_max;

Rejected(:,10) = P./C < config.P_to_C_min;
Rejected(:,11) = P./C > config.P_to_C_max;

Rejected(:,12) = O./P < config.O_to_P_min;
Rejected(:,13) = O./P > config.O_to_P_max;

Rejected(:,14) = DBE < config.DBE_min;
Rejected(:,15) = DBE > config.DBE_max;

Rejected(:,16) = AImod > 1;

Rejected(:,17) = K+Na >1; % Removes any formulas with both K and Na assigned

if strcmpi(config.Ionization_Mode,'negative')
    if config.Heteroelement_Halogen
        Rejected(:,18) = (O < 1 & S < 1 & E < 1); % Removes any formulas with O=0 & S=0 & E=0(i.e., CH, CHN, CHP)
    else
        Rejected(:,18) = (O < 1 & S < 1 ); % Removes any formulas with O=0 & S=0 (i.e., CH, CHN, CHP, CHE)
    end
elseif strcmpi(config.Ionization_Mode,'positive')
    Rejected(:,18) = (O < 1 & S < 1 & N < 0); % Removes any formulas with O=0 & S=0 & N=0(i.e., CH, CHP, CHE)
else
    error('Error! Need to select the correct ionization mode (positive or negative) in the configuration file!')
end
    
Data_Stage2= [Index mz Magnitude SN EstimC C H O N S P E ExactMass e OC HC DBE DBEC AImod];

% Define new format (internal to the code!)
global format 

format.Column_Index=1;
format.Column_mz=2;
format.Column_Magnitude=3;
format.Column_SN=4;
format.Column_EstimC=5;
format.Column_C=6;
format.Column_H=7;
format.Column_O=8;
format.Column_N=9;
format.Column_S=10;
format.Column_P=11;
format.Column_E=12;
format.Column_ExactMass=13;
format.Column_error=14;
format.Column_OC=15;
format.Column_HC=16;
format.Column_DBE=17;
format.Column_DBEC=18;
format.Column_AImod=19;
format.elements = ["C", "H", "O", "N", "S", "P", "E"];

% Split the imported data (Data_Stage2) into two fractions:
Data_Stage2_Rejected = Data_Stage2(any(Rejected')',:); 
Data_Stage2_Refined = Data_Stage2(~any(Rejected')',:);

titles={'m/z' 'Magnitude' 'S/N' 'Estim C#' 'Type' 'C' 'H' 'O' 'N' 'S' 'P' 'E' 'ExactMass' 'Assignment Error (ppm)' 'O/C' 'H/C' 'DBE' 'DBE/C' 'AImod'};

% Export
if ~isempty(Data_Stage2_Rejected) % There were some rejected
    xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Rejected using Elem.Constraints','A1');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage2_Rejected),'Rejected using Elem.Constraints','A2');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],{'Formulas eliminated using the elemental constraints (Stubbins et al. 2010 and others). H is corrected!'},'Rejected using Elem.Constraints','U2');
    
    xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Refined using Elem.Constraints','A1');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage2_Refined),'Refined using Elem.Constraints','A2');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],{'Data that fits the elemental constraints (Stubbins et al. 2010 and others), H is corrected, ready for further refienement!'},'Refined using Elem.Constraints','U2');
else % None were rejected - keep the original list
    xlswrite([filename(1:end-6) '_Processing.xlsx'],{'No formulas were rejected using the elemental constraints (Stubbins et al. 2010 and others)'},'Rejected using Elem.Constraints','A1');
    
    xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Refined using Elem.Constraints','A1');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage2_Refined),'Refined using Elem.Constraints','A2');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],{'Data that fits the elemental constraints (Stubbins et al. 2010 and others), H is corrected, ready for further refienement!'},'Refined using Elem.Constraints','U2');
end

%% Stage 3: Formula Refinement using Isotopic Filter and selection of Unique Formulas

clearvars -except filename config format titles RefinementData RefinementQC Data_Stage2_Refined 

Data_Stage3 = sortrows(Data_Stage2_Refined, format.Column_mz);

% Isotopic filter
Data_Stage3_Isotopes = [];
if config.Filter_Isotopes 
    EstimC = Data_Stage3(:,format.Column_mz:format.Column_EstimC);
    
    EstimC = EstimC(EstimC(:,3) > config.Filter_Isotopes_SNthreshold,[1,2,3,4]);    % Pick m/z above the S/N threshold
    EstimC = EstimC(EstimC(:,4) ~= 0,[1,4,4]);                                      % Pick m/z containing a C-estimate
    EstimC(:,2) = EstimC(:,2) - config.Filter_Isotopes_Range;
    EstimC(:,3) = EstimC(:,3) + config.Filter_Isotopes_Range;
    EstimC(:,2:3) = floor(EstimC(:,2:3));
    
    i = 1; j = 1;
    while i <= size(EstimC,1) && j <= size(Data_Stage3, 1)
        mz_BeforeAssignment = EstimC(i, 1);
        mz_AfterAssignment = Data_Stage3(j, format.Column_mz);
        if intersect(mz_BeforeAssignment,mz_AfterAssignment) > 0
            AssignedCnumber = Data_Stage3(j, format.Column_C);
            CalcCnumber_lower = EstimC(i, 2);
            CalcCnumber_higher = EstimC(i, 3);
            if (AssignedCnumber >= CalcCnumber_lower && AssignedCnumber <= CalcCnumber_higher)
                Data_Stage3_Isotopes = [Data_Stage3_Isotopes ; Data_Stage3(j,:)];
            end
            j = j+1;
        elseif mz_BeforeAssignment < mz_AfterAssignment
            i = i+1;
        elseif mz_BeforeAssignment > mz_AfterAssignment
            j = j+1;
        end
    end
    
    if ~isempty(Data_Stage3_Isotopes)
        xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Isotopic Filter','A1');
        xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage3_Isotopes),'Isotopic Filter','A2');
    else
        xlswrite([filename(1:end-6) '_Processing.xlsx'],{'No formulas were refined with the isotopic filter'},'Isotopic Filter','A1');
    end
end

% Uniqueness filter
Previous_Index = 0;
Data_Stage3_Unique = [];
if config.Filter_Unique
    for i=1:size(Data_Stage3,format.Column_Index)
        if i < size(Data_Stage3,format.Column_Index) - 1
            if Data_Stage3(i,format.Column_Index) ~= Data_Stage3(i+1,format.Column_Index) && Data_Stage3(i,format.Column_Index) ~= Previous_Index
                Data_Stage3_Unique = [Data_Stage3_Unique; Data_Stage3(i,:)];
            end
        elseif Data_Stage3(i,format.Column_Index) ~= Previous_Index
            Data_Stage3_Unique = [Data_Stage3_Unique; Data_Stage3(i,:)];
        end
        Previous_Index = Data_Stage3(i,format.Column_Index);
    end
end

if ~config.Filter_Unique && ~config.Filter_Isotopes
    error('Need to enable EITHER "Accepting all unique formulas" OR "Isotopic selection" in  config file!')
end

if ~isempty(Data_Stage3_Unique)
    xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Unique','A1');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage3_Unique),'Unique','A2');
else
    xlswrite([filename(1:end-6) '_Processing.xlsx'],{'No unique formulas found'},'Unique','A1');
end

if isempty(Data_Stage3_Isotopes) && isempty(Data_Stage3_Unique)
    Data_Stage3_Refined=Data_Stage3;
else
    Data_Stage3_Refined = [Data_Stage3_Isotopes ; Data_Stage3_Unique];
    % Remove duplicates - if any unique formulas happen to be ALSO selected by the isotopic filter, 
    % duplicate entries will exist at this point. Here we remove them.
    Data_Stage3_Refined = unique(Data_Stage3_Refined, 'rows');
end

%% Stage 4: Formula refinement using KMD Filter

clearvars -except filename config format titles RefinementData RefinementQC Data_Stage3 Data_Stage3_Refined    

Data_Stage4 = sortrows(Data_Stage3, format.Column_ExactMass);                 % All data
Data_Stage4_Refined = sortrows(Data_Stage3_Refined, format.Column_ExactMass); % Formulas of high confidence

if config.Filter_KMD
    % remove formulas of high confidence from all formulas -> identifies which formulas can be refined
    indices_to_remove = [];
    i = 1; j = 1;
    
    % Identify the indices of "formulas of high confidence" in "all data"
    while i <= size(Data_Stage4_Refined,1) && j <= size(Data_Stage4, 1) 
        if intersect(Data_Stage4_Refined(i,format.Column_ExactMass),Data_Stage4(j,format.Column_ExactMass)) > 0
            indices_to_remove = [indices_to_remove ; j];
            j = j+1;
        elseif Data_Stage4_Refined(i,format.Column_ExactMass) < Data_Stage4(j,format.Column_ExactMass)
            i = i+1;
        else
            j = j+1;
        end
    end
    
    Data_Stage4_Refinable = Data_Stage4(setdiff(1:size(Data_Stage4,1),indices_to_remove),:);
    search_values = getCompositionAndKMD(Data_Stage4_Refinable,config);
    kmd_counts = getKMDCounts(Data_Stage4_Refined,config);

    should_refine = true;
    while should_refine
        tic
        new_candidate_ids = kmd_counts.FindMatchingCandidates(search_values, config);
        removed_candidates = search_values.RemoveCandidates(new_candidate_ids);
        kmd_counts.UpdateKMDValues(removed_candidates, config);
        new_candidates = Data_Stage4_Refinable(new_candidate_ids,:);
        Data_Stage4_Refined = [Data_Stage4_Refined ; new_candidates];

        should_refine = size(new_candidate_ids,1) > 0;
    end

    % Export all formulas of high confidence (isotope and/or unique + KMD refined)
    if ~isempty(Data_Stage4_Refined)
        xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Refined w KMD','A1');
        xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage4_Refined),'Refined w KMD','A2');
    else
        xlswrite([filename(1:end-6) '_Processing.xlsx'],{'No formulas were refined with KMD'},'Refined w KMD','A1');
    end
    
else
    Data_Stage4_Refinable=Data_Stage4;
end

Data_Stage4_Refined   = sortrows(Data_Stage4_Refined,format.Column_Index);
Data_Stage4_Refinable = sortrows(Data_Stage4_Refinable,format.Column_Index);

% Identify peaks for which we have not yet identified a good candidate molecule (refinable) 
% Also identify peaks which we have already identified good candidates for (and discard the other ambiguous assignments)
i = 1; j = 1; idx = [];
final_size = size(Data_Stage4_Refined, 1);
while i <= final_size && j <= size(Data_Stage4_Refinable, 1)
    index_refined = Data_Stage4_Refined(i, format.Column_Index);
    index_refinable = Data_Stage4_Refinable(j, format.Column_Index);
    if index_refined == index_refinable
        j = j+1;
    elseif index_refined < index_refinable
        i = i+1;
    elseif index_refined > index_refinable
        idx = [idx ; j];
        Data_Stage4_Refined = [Data_Stage4_Refined ; Data_Stage4_Refinable(j,:)];
        j = j+1;
    end
end

% Removes any ambiguous formulas to those found in the refined list(Data_Stage4_Refined).

% If KMD filter is disabled, it will only list ambiguous formulas to those selected by the isotope filter
% If KMD filter is enabled, it will list ambiguous formulas to those selected by the isotope and KMD filters
Data_Stage4_Rejected = Data_Stage4_Refinable(setdiff(1:size(Data_Stage4_Refinable,1),idx),:);

if ~isempty(Data_Stage4_Rejected)
    xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Rejected at Middle Stage','A1');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage4_Rejected),'Rejected at Middle Stage','A2');
else
    xlswrite([filename(1:end-6) '_Processing.xlsx'],{'No formulas were rejected with the isotope filter and/or KMD'},'Rejected at Middle Stage','A1');
end


xlswrite([filename(1:end-6) '_Processing.xlsx'],titles, 'Middle Stage','A1');
xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage4_Refined), 'Middle Stage','A2');

%% Stage 5: Formula refinement using Compositional Filter

clearvars -except filename config format titles RefinementData RefinementQC Data_Stage4_Refined 

Data_Stage5=sortrows(Data_Stage4_Refined,format.Column_Index); 

if config.Filter_Composition
    simple_batch = [];
    current_id = -1;
    indices = [];
    for i=1:size(Data_Stage5,1)
        if Data_Stage5(i,1) ~= current_id
            if size(simple_batch,1) > 0
                simple = [];
                everything = [];
                for j=1:size(simple_batch,1)
                    everything = [everything ; simple_batch{j, 1}];
                    if simple_batch{j, 2}
                        simple = [simple ; simple_batch{j, 1}];
                    end
                end
                if size(simple,1) > 0
                    indices = [indices ; simple];
                else
                    indices = [indices ; everything];
                end
            end
            simple_batch = [];
            current_id = Data_Stage5(i,1);
        end
        simple_batch = [simple_batch ; {i, evaluateSimple(Data_Stage5(i,:))}];
    end
    
    simple = [];
    everything = [];
    for i=1:size(simple_batch,1)
        everything = [everything ; simple_batch{i, 1}];
        if simple_batch{i, 2}
            simple = [simple ; simple_batch{i, 1}];
        end
    end
    
    if size(simple,1) > 0
        indices = [indices ; simple];
    else
        indices = [indices ; everything];
    end
    
    Data_Stage5_Rejected = Data_Stage5(setdiff(1:size(Data_Stage5,1),indices),:);
    Data_Stage5_Refined = Data_Stage5(indices,:);

    if ~isempty(Data_Stage5_Rejected)
        xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Rejected w Compositional Filter','A1');
        xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage5_Rejected),'Rejected w Compositional Filter','A2');

        xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Refined w Compositional Filter','A1');
        xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage5_Refined),'Refined w Compositional Filter','A2');      
    end
else
    Data_Stage5_Refined=Data_Stage5; % If composition filter is disabled, this transfers the original data into the next stage
end
%% Stage 6: Formula refinement using Error Filter

clearvars -except filename config format titles RefinementData RefinementQC Data_Stage5_Refined 

Data_Stage6=sortrows(Data_Stage5_Refined,format.Column_Index); 

if config.Filter_Error
    error_batch = [];
    current_id = -1;
    indices = [];
    for i=1:size(Data_Stage6,format.Column_Index)
        if Data_Stage6(i,format.Column_Index) ~= current_id
            if size(error_batch,1) > 0
                error_batch = sortrows(error_batch, 2); % Sort so lowest is on top (1,2)
                best_error = error_batch{1,2};
                j = 1;
                while j <= size(error_batch,1) && ~isempty(intersect(error_batch{j,2},best_error))
                    indices = [indices ; error_batch{j,1}];
                    j = j+1;
                end
            end
            error_batch = [];
            current_id = Data_Stage6(i,format.Column_Index);
        end
        error_batch = [error_batch ; {i, abs(Data_Stage6(i,format.Column_error))}];
    end
    error_batch = sortrows(error_batch, 2);
    best_error = error_batch{1,2};
    j = 1;
    while j <= size(error_batch,1) && ~isempty(intersect(error_batch{j,2},best_error))
        indices = [indices ; error_batch{j,1}];
        j = j+1;
    end
    Data_Stage6_Rejected = Data_Stage6(setdiff(1:size(Data_Stage6,1),indices),:);
    Data_Stage6_Refined = Data_Stage6(indices,:);

    if ~isempty(Data_Stage6_Rejected)
        xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Rejected w Error Filter','A1');
        xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage6_Rejected),'Rejected w Error Filter','A2');
    end
    
    xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Refined w Error Filter','A1');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage6_Refined),'Refined w Error Filter','A2');
else
    Data_Stage6_Refined=Data_Stage6; % If the error filter is disabled, this transfers the original data into the next stage
end
%% Stage 7: Final Refinement

clearvars -except filename config format titles RefinementData RefinementQC Data_Stage6_Refined  

Data_Stage7=sortrows(Data_Stage6_Refined,format.Column_Index);
Data_Stage7_Refined1=[];
Data_Stage7_Rejected1=[];

% Final refinement 1: Remove any residual duplicates, two or more identical indexes (i.e., two identical m/z values)
for i=1:size(Data_Stage7,1)
    [r,~]=find(Data_Stage7(:,format.Column_Index)==Data_Stage7(i,format.Column_Index)); 
    if length(r) == 1
        Data_Stage7_Refined1=[Data_Stage7_Refined1;Data_Stage7(i,:)];
    else
        Data_Stage7_Rejected1=[Data_Stage7_Rejected1;Data_Stage7(i,:)];
    end
end

Data_Stage7_Refined2=[];
Data_Stage7_Refinable2=[];
Data_Stage7_Rejected2=[];

% Final refinement 2: Remove exact mass duplicates
for i=1:size(Data_Stage7_Refined1,1)
    [r,~]=find(round(Data_Stage7_Refined1(:,format.Column_ExactMass),config.Precision) == round(Data_Stage7_Refined1(i,format.Column_ExactMass),config.Precision));
    if length(r) == 1       % Unambiguous
        Data_Stage7_Refined2=[Data_Stage7_Refined2;Data_Stage7_Refined1(i,:)];
    elseif length(r) == 2   % Ambiguous, 2 options
        Data_Stage7_Refinable2=[Data_Stage7_Refinable2;Data_Stage7_Refined1(i,:)];
    else
        Data_Stage7_Rejected2=[Data_Stage7_Rejected2;Data_Stage7_Refined1(i,:)];
        disp(['ExactMass = ' char(string(Data_Stage7_Refined1(i,format.Column_ExactMass))) ' is assigned to more than two peaks - inspect spectrum!'])
    end
end

Data_Stage7_Refined3=[];
Data_Stage7_Rejected3=[];

% Pick one of the ExactMass Duplicates
if ~isempty(Data_Stage7_Refinable2)
    for i=1:size(Data_Stage7_Refinable2,1)
        [r,~]=find(round(Data_Stage7_Refinable2(:,format.Column_ExactMass),config.Precision) == round(Data_Stage7_Refinable2(i,format.Column_ExactMass),config.Precision)); 

        % In case of shoulder peaks (one of the peaks is bigger with more than 20% & bigger peak has a lower error):
        % DO: Accept the larger peak  
        if Data_Stage7_Refinable2(r(1),format.Column_Magnitude) > Data_Stage7_Refinable2(r(2),format.Column_Magnitude) && ...
                ((Data_Stage7_Refinable2(r(1),format.Column_Magnitude))-(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)))/(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)) >= 0.20 ... 
                && abs(Data_Stage7_Refinable2(r(1),format.Column_error)) < abs(Data_Stage7_Refinable2(r(2),format.Column_error)) ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_mz)-Data_Stage7_Refinable2(r(2),format.Column_mz)) < 1
            
            Data_Stage7_Refined3=[Data_Stage7_Refined3;Data_Stage7_Refinable2(r(1),:)];
            Data_Stage7_Rejected3=[Data_Stage7_Rejected3;Data_Stage7_Refinable2(r(2),:)];
            
        elseif Data_Stage7_Refinable2(r(1),format.Column_Magnitude) < Data_Stage7_Refinable2(r(2),format.Column_Magnitude) && ...
                ((Data_Stage7_Refinable2(r(1),format.Column_Magnitude))-(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)))/(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)) >= 0.20 ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_error)) > abs(Data_Stage7_Refinable2(r(2),format.Column_error)) ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_mz)-Data_Stage7_Refinable2(r(2),format.Column_mz)) < 1
            
            Data_Stage7_Refined3=[Data_Stage7_Refined3;Data_Stage7_Refinable2(r(2),:)];
            Data_Stage7_Rejected3=[Data_Stage7_Rejected3;Data_Stage7_Refinable2(r(1),:)];

        % In case of split peaks (difference in spectral magnitude is less than 20%):
        % DO: Accept the peak with smaller error       
        elseif ((Data_Stage7_Refinable2(r(1),format.Column_Magnitude))-(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)))/(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)) < 0.20 ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_error)) < abs(Data_Stage7_Refinable2(r(2),format.Column_error)) ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_mz)-Data_Stage7_Refinable2(r(2),format.Column_mz)) < 1
            
            Data_Stage7_Refined3=[Data_Stage7_Refined3;Data_Stage7_Refinable2(r(1),:)];
            Data_Stage7_Rejected3=[Data_Stage7_Rejected3;Data_Stage7_Refinable2(r(2),:)];
            
        elseif ((Data_Stage7_Refinable2(r(1),format.Column_Magnitude))-(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)))/(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)) < 0.20 ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_error)) > abs(Data_Stage7_Refinable2(r(2),format.Column_error)) ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_mz)-Data_Stage7_Refinable2(r(2),format.Column_mz)) < 1
            
            Data_Stage7_Refined3=[Data_Stage7_Refined3;Data_Stage7_Refinable2(r(2),:)];
            Data_Stage7_Rejected3=[Data_Stage7_Rejected3;Data_Stage7_Refinable2(r(1),:)];

        % In case of an intense spike peaks (one of the peaks is bigger with more than 20% & bigger peak has a bigger error)
        % DO: Discard both and output a message
        elseif Data_Stage7_Refinable2(r(1),format.Column_Magnitude) > Data_Stage7_Refinable2(r(2),format.Column_Magnitude) ...
                && ((Data_Stage7_Refinable2(r(1),format.Column_Magnitude))-(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)))/(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)) >= 0.20 ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_error)) > abs(Data_Stage7_Refinable2(r(2),format.Column_error)) ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_mz)-Data_Stage7_Refinable2(r(2),format.Column_mz)) < 1
            
            Data_Stage7_Refined3=Data_Stage7_Refined3;
            Data_Stage7_Rejected3=[Data_Stage7_Rejected3;Data_Stage7_Refinable2(r(1),:);Data_Stage7_Refinable2(r(2),:)];
            disp(['Peaks with ExactMass = ' char(string(Data_Stage7_Refined1(r(1),format.Column_ExactMass))) ' need inspection and manual assignment!'])
            
        elseif Data_Stage7_Refinable2(r(1),format.Column_Magnitude) < Data_Stage7_Refinable2(r(2),format.Column_Magnitude) ...
                && ((Data_Stage7_Refinable2(r(1),format.Column_Magnitude))-(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)))/(Data_Stage7_Refinable2(r(2),format.Column_Magnitude)) >= 0.20 ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_error)) < abs(Data_Stage7_Refinable2(r(2),format.Column_error)) ...
                && abs(Data_Stage7_Refinable2(r(1),format.Column_mz)-Data_Stage7_Refinable2(r(2),format.Column_mz)) < 1
            
            Data_Stage7_Refined3=Data_Stage7_Refined3;
            Data_Stage7_Rejected3=[Data_Stage7_Rejected3;Data_Stage7_Refinable2(r(1),:);Data_Stage7_Refinable2(r(2),:)];
            disp(['Peaks with ExactMass = ' char(string(Data_Stage7_Refined1(r(2),format.Column_ExactMass))) ' need inspection and manual assignment!'])
        
        % In case of two different ion species were detected: 
        % [M+H]+ and [M+Na]+; [M+H]+ and [M+K]+; [M+K]+ and [M+Na]+; or other combos
        % DO: Choose the smaller option (e.g., if [M+H]+ and [M+Na]+ are
        % detected, pick [M+H]+ 
        elseif abs(Data_Stage7_Refinable2(r(1),format.Column_mz)-Data_Stage7_Refinable2(r(2),format.Column_mz)) > 10 ...
                && Data_Stage7_Refinable2(r(1),format.Column_mz) < Data_Stage7_Refinable2(r(2),format.Column_mz) 
            
            Data_Stage7_Refined3=[Data_Stage7_Refined3;Data_Stage7_Refinable2(r(1),:)];
            Data_Stage7_Rejected3=[Data_Stage7_Rejected3;Data_Stage7_Refinable2(r(2),:)];
            
        elseif abs(Data_Stage7_Refinable2(r(1),format.Column_mz)-Data_Stage7_Refinable2(r(2),format.Column_mz)) > 10 ...
                && Data_Stage7_Refinable2(r(1),format.Column_mz) > Data_Stage7_Refinable2(r(2),format.Column_mz) 
        
            Data_Stage7_Refined3=[Data_Stage7_Refined3;Data_Stage7_Refinable2(r(2),:)];
            Data_Stage7_Rejected3=[Data_Stage7_Rejected3;Data_Stage7_Refinable2(r(1),:)];        
        
        else
            ;
        end
    end
end

Data_Stage7_Refined=[Data_Stage7_Refined2;Data_Stage7_Refined3];
Data_Stage7_Rejected=[Data_Stage7_Rejected1;Data_Stage7_Rejected2;Data_Stage7_Rejected3];

if ~isempty(Data_Stage7_Rejected)
    xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Rejected w Final Refinement','A1');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage7_Rejected),'Rejected w Final Refinement','A2');
else
    xlswrite([filename(1:end-6) '_Processing.xlsx'],{'No furter refinement was needed'},'Rejected w Final Refinement','A1');
end
 
% Identify formulas that were not assigned
Assigned=intersect(RefinementData(:,1),Data_Stage7_Refined(:,format.Column_mz));
[~,~,index_c]=intersect(Assigned,RefinementData(:,1));
index_u=[1:1:size(RefinementData(:,1),1)]'; index_u(index_c)=[];
Unassigned=RefinementData(index_u,:);

if ~isempty(Unassigned)
    xlswrite([filename(1:end-6) '_Processing.xlsx'],{'m/z' 'Magnitude','S/N','Estim C#'},'Not Assigned','A1');
    xlswrite([filename(1:end-6) '_Processing.xlsx'],Unassigned,'Not Assigned','A2');   
else
    xlswrite([filename(1:end-6) '_Processing.xlsx'],{'No formulas were not assignable (assignment = 100%)'},'Rejected w Final Refinement','A1');
end

xlswrite([filename(1:end-6) '_Processing.xlsx'],titles,'Final Refinement','A1');
xlswrite([filename(1:end-6) '_Processing.xlsx'],FormatDataForExport(Data_Stage7_Refined),'Final Refinement','A2');

xlswrite([filename(1:end-17) '_Final.xlsx'],titles,'Sheet1','A1');
xlswrite([filename(1:end-17) '_Final.xlsx'],FormatDataForExport(Data_Stage7_Refined),'Sheet1','A2');
    
%% Stage 8: Quality Control
clearvars -except filename config format titles RefinementData RefinementQC Data_Stage7_Refined Unassigned

figure 
set(gcf,'color',[0.85 0.85 0.85]);
set(gcf,'WindowState','maximized')
    subplot(2,3,[1,4]) %vK evaluating S/N
        hold on
        title('Van Krevelen Diagram');
        xlabel('O/C'); ylabel('H/C');
        xlim([0 1.2]); ylim([0 2.5]);
        figure_matrix=[Data_Stage7_Refined(:,format.Column_OC),Data_Stage7_Refined(:,format.Column_HC),Data_Stage7_Refined(:,format.Column_SN)];
        figure_matrix=sortrows(figure_matrix,3,'descend');
        figure_matrix=figure_matrix(10:end,:); % Remove top 10 to avoid outliers skewing the color scheme
        scatter3(figure_matrix(:,1),figure_matrix(:,2),figure_matrix(:,3),70,figure_matrix(:,3),'filled','o');
        colormap(jet);
        handle_colorbar=colorbar;
        set(get(handle_colorbar,'title'),'string','S/N','fontsize',13);
        view(0,90)     
        h1=refline(-1,2);           % AImod=0
        h2=refline(-0.3845,1.0584); % AImod=0.5
        h3=refline(-0.3740,0.7551); % AImod=0.67
        h1.Color='k';h2.Color='k';h3.Color='k';
        h1.LineWidth=1.5;h2.LineWidth=1.5;h3.LineWidth=1.5;
        set(gca,'fontweight','bold')
        hold off
    subplot(2,3,2) % Histogram for S/N 
        hold on
        title('Histogram S/N');
        xlabel('S/N'); ylabel('Number of Formulas');
        histfit(Data_Stage7_Refined(:,format.Column_SN),150);
        dist1=fitdist(Data_Stage7_Refined(:,format.Column_SN),'Normal');
        if dist1.sigma ~=0
            title(['Histogram Signal-to-Noise (' num2str(dist1.mu) ' ± ' num2str(dist1.sigma) ')']);
        else
            title('Histogram Signal-to-Noise (N/A)');
        end
        xlim([-50 300])
        set(gca,'fontweight','bold')
        hold off
    subplot(2,3,5) % Histogram for Assignment Errors
        hold on
        xlabel('Assignment Error (ppm)'); ylabel('Number of Formulas');
        histfit(Data_Stage7_Refined(:,format.Column_error),21);
        dist2=fitdist(Data_Stage7_Refined(:,format.Column_error),'Normal');
        title(['Histogram Assignment Error (' num2str(round(dist2.mu,4)) ' ± ' num2str(round(dist2.sigma,4)) ' ppm) ']);
        set(gca,'fontweight','bold')
        hold off
    subplot(2,3,3) % Evaluation of C-numbers
        hold on
        xlabel('Estimated C-number from ^{13}C/^{12}C ratio'); ylabel('Assigned C-number');
        Cnumbers=[Data_Stage7_Refined(:,format.Column_EstimC),Data_Stage7_Refined(:,format.Column_C)];
        Cnumbers(Cnumbers == 0) = NaN;
        Cnumbers(any(isnan(Cnumbers), 2), :) = [];
        plot(Cnumbers(:,1),Cnumbers(:,2),'Marker','.','MarkerSize',13,'Color',[0.0745098039215686 0.623529411764706 1],'LineStyle','none');
        linefit=polyfit(Cnumbers(:,1),Cnumbers(:,2),1);
        R=corrcoef(Cnumbers(:,1),Cnumbers(:,2));
        Rsq = R(1,2).^2;
        plot([min([xlim ylim]) max([xlim ylim])], [min([xlim ylim]) max([xlim ylim])], '-r','LineWidth',3); % 1:1 line
        h4=refline(linefit(1),linefit(2));
        set(h4,'LineWidth',3,'Color','k','LineStyle','--') % Linear Fit
        title(['C Assessment (m = ' num2str(round(linefit(1),3)) ', b = ' num2str(round(linefit(2),3)) ', R^2 = ' num2str(round(Rsq,3)) ')']);
        legend({'C numbers','One-to-One Line','Linear Fit'},'Location','southeast')
        set(gca,'fontweight','bold')
        hold off
    subplot(2,3,6) % Table with stats for the processing
        hold on
        ylim([0 1.4]); xlim([0 1]);
        set(gca,'YTickLabel',[],'XTickLabel',[],'XColor','none','YColor','none');
        text(0.01,1.25,'Number of Peaks','FontWeight','bold','FontSize',12)
        text(0.01,1.15,'Blank Peaks','FontWeight','bold','FontSize',12)
        text(0.01,1.05,'Salt Peaks','FontWeight','bold','FontSize',12)
        text(0.01,0.95,'Doubly Charged Peaks','FontWeight','bold','FontSize',12)
        text(0.01,0.85,'^{13}C Peaks','FontWeight','bold','FontSize',12)
        text(0.01,0.75,'^{34}S Peaks','FontWeight','bold','FontSize',12)
        text(0.01,0.65,'^{54}Fe Peaks','FontWeight','bold','FontSize',12)
        text(0.01,0.55,'^{37}Cl Peaks','FontWeight','bold','FontSize',12)
        text(0.01,0.45,'^{200}Hg Peaks','FontWeight','bold','FontSize',12)
        text(0.01,0.35,'Rejected Peaks','FontWeight','bold','FontSize',12)
        text(0.01,0.25,'Refined Peaks','FontWeight','bold','FontSize',12)
        text(0.01,0.15,'Unassignable Peaks','FontWeight','bold','FontSize',12)
        
        text(0.4646,1.35,'# of peaks   % Num  % Magn','FontWeight','bold','FontSize',12)
        annotation('line',[0.7876 0.7876],[0.4255 0.12],'LineWidth',2);
        annotation('line',[0.8362 0.8362],[0.4255 0.12],'LineWidth',2);
        annotation('line',[0.8709 0.8709],[0.4255 0.12],'LineWidth',2);
        annotation('line',[0.7876 0.9035],[0.4255 0.4255],'LineWidth',2);
        
        text(0.5173,1.25,num2str(RefinementQC(1,1)),'FontSize',12); % Number of formulas
        text(0.5173,1.15,num2str(RefinementQC(2,1)),'FontSize',12);
        text(0.5173,1.05,num2str(RefinementQC(3,1)),'FontSize',12);
        text(0.5173,0.95,num2str(RefinementQC(4,1)),'FontSize',12);
        text(0.5173,0.85,num2str(RefinementQC(5,1)),'FontSize',12);
        text(0.5173,0.75,num2str(RefinementQC(6,1)),'FontSize',12);
        text(0.5173,0.65,num2str(RefinementQC(7,1)),'FontSize',12);
        text(0.5173,0.55,num2str(RefinementQC(8,1)),'FontSize',12);
        text(0.5173,0.45,num2str(RefinementQC(9,1)),'FontSize',12);
        text(0.5173,0.35,num2str(RefinementQC(10,1)),'FontSize',12);
        text(0.5173,0.25,num2str(RefinementQC(11,1)),'FontSize',12);
        text(0.5173,0.15,num2str(size(Unassigned,1)),'FontSize',12);
        
        text(0.7124,1.25,num2str(round(RefinementQC(1,2),2)),'FontSize',12); % Percentages, number-based
        text(0.7124,1.15,num2str(round(RefinementQC(2,2),2)),'FontSize',12);
        text(0.7124,1.05,num2str(round(RefinementQC(3,2),2)),'FontSize',12);
        text(0.7124,0.95,num2str(round(RefinementQC(4,2),2)),'FontSize',12);
        text(0.7124,0.85,num2str(round(RefinementQC(5,2),2)),'FontSize',12);
        text(0.7124,0.75,num2str(round(RefinementQC(6,2),2)),'FontSize',12);
        text(0.7124,0.65,num2str(round(RefinementQC(7,2),2)),'FontSize',12);
        text(0.7124,0.55,num2str(round(RefinementQC(8,2),2)),'FontSize',12);
        text(0.7124,0.45,num2str(round(RefinementQC(9,2),2)),'FontSize',12);
        text(0.7124,0.35,num2str(round(RefinementQC(10,2),2)),'FontSize',12);
        text(0.7124,0.25,num2str(round(RefinementQC(11,2),2)),'FontSize',12);
        text(0.7124,0.15,num2str(round(size(Unassigned,1)*100/RefinementQC(11,1),2)),'FontSize',12);
        
        text(0.8831,1.25,num2str(round(RefinementQC(1,3),2)),'FontSize',12); % Percentages, magnitude-based
        text(0.8831,1.15,num2str(round(RefinementQC(2,3),2)),'FontSize',12);
        text(0.8831,1.05,num2str(round(RefinementQC(3,3),2)),'FontSize',12);
        text(0.8831,0.95,num2str(round(RefinementQC(4,3),2)),'FontSize',12);
        text(0.8831,0.85,num2str(round(RefinementQC(5,3),2)),'FontSize',12);
        text(0.8831,0.75,num2str(round(RefinementQC(6,3),2)),'FontSize',12);
        text(0.8831,0.65,num2str(round(RefinementQC(7,3),2)),'FontSize',12);
        text(0.8831,0.55,num2str(round(RefinementQC(8,3),2)),'FontSize',12);
        text(0.8831,0.45,num2str(round(RefinementQC(9,3),2)),'FontSize',12);
        text(0.8831,0.35,num2str(round(RefinementQC(10,3),2)),'FontSize',12);
        text(0.8831,0.25,num2str(round(RefinementQC(11,3),2)),'FontSize',12);
        text(0.8831,0.15,num2str(round(sum(Unassigned(:,2))*100/sum(RefinementData(:,2)),2)),'FontSize',12);

        % Assigned peaks computation
        if round(size(Data_Stage7_Refined,1)*100/RefinementQC(11,1),2) > 80 || round(sum(Data_Stage7_Refined(:,format.Column_Magnitude))*100/sum(RefinementData(:,2)),2) > 80
            text(0.01,0.05,'Assigned Peaks','FontWeight','bold','FontSize',12,'Color',[0.25,0.67,0.03])
            text(0.5173,0.05,num2str(size(Data_Stage7_Refined,1)),'FontWeight','bold','FontSize',12,'Color',[0.25,0.67,0.03]);
            text(0.7124,0.05,num2str(round(size(Data_Stage7_Refined,1)*100/RefinementQC(11,1),2)),'FontWeight','bold','FontSize',12,'Color',[0.25,0.67,0.03]);
            text(0.8831,0.05,num2str(round(sum(Data_Stage7_Refined(:,format.Column_Magnitude))*100/sum(RefinementData(:,2)),2)),'FontWeight','bold','FontSize',12,'Color',[0.25,0.67,0.03]);
        else
            text(0.01,0.05,'Assigned Peaks','FontWeight','bold','FontSize',12,'Color','r')
            text(0.5173,0.05,num2str(size(Data_Stage7_Refined,1)),'FontWeight','bold','FontSize',12,'Color','r');
            text(0.7124,0.05,num2str(round(size(Data_Stage7_Refined,1)*100/RefinementQC(11,1),2)),'FontWeight','bold','FontSize',12,'Color','r');
            text(0.8831,0.05,num2str(round(sum(Data_Stage7_Refined(:,format.Column_Magnitude))*100/sum(RefinementData(:,2)),2)),'FontWeight','bold','FontSize',12,'Color','r');
        end
print(gcf,['FTMS Processing_' filename(1:end-17)],'-dpng','-r300');
disp(['Finished refining the formulas of the peak list ' char(filename) ' (' num2str(toc) ' seconds)'])
end

%% Internal functions

% Function for Stage1: gets numerical info from header lines beginning with "Peak Report" 
function [hdrdat]=retheader(str) 
fstring='Peak Report for m/z = %f +/-%f ppm     Peak Height = %f';  % define string format to locate numerical values
hdrdat=sscanf(str,fstring);  % pull numerical values from header line
end

% Function for Stage1: gets numerical data from lined containing elemental values (C,H,O,N,S,P,E,K,Na)
function [dat]=retdata(str) 
% gets numerical values from data string - this line must be changed if elemental output is changed ie elements added, removed, or order changed
%          C         H         N          O          S           P         Cl         K          Na          m/z        error   
datstr={str(2:5),str(7:10),str(12:15),str(17:20),str(22:25),str(27:30),str(38:40),str(42:45),str(48:50),str(57:66),str(68:71)}; 
dat=str2double(datstr);  % convert data from string to number
end

% Function for Stage4: Counts the number of points on a KMD series
function kmd_totals = getKMDCounts(compound_list,config)
    composition_kmd = getCompositionAndKMD(compound_list,config);
    kmd_totals = FTMS_KMD_Collection(composition_kmd, config);
end

% Function for Stage4: Gets the string composition and kmd values of a given matrix
function composition_kmd = getCompositionAndKMD(compound_list,config)
    composition_kmd = FTMS_CandidateCollection(compound_list,config);
end

% Function for Stage 5: Evaluation of formula composition (simple vs complex)
function is_simple = evaluateSimple(compound)
    global format
    persistent C_index H_index O_index N_index P_index S_index other_indices;
    if isempty(C_index)
        C_index = format.Column_C - 1 + find(format.elements=="C");
        H_index = format.Column_C - 1 + find(format.elements=="H");
        O_index = format.Column_C - 1 + find(format.elements=="O");
        N_index = format.Column_C - 1 + find(format.elements=="N");
        P_index = format.Column_C - 1 + find(format.elements=="P");
        S_index = format.Column_C - 1 + find(format.elements=="S");
        chonps_indices = [C_index H_index O_index N_index P_index S_index];
        other_indices = setdiff(format.Column_C:(format.Column_C + size(format.elements,2) - 1),chonps_indices);
    end
    sum_others = sum(compound(1,other_indices));
    if sum_others ~= 0
        is_simple = false;
    else
        if compound(1,C_index) == 0 || compound(1,H_index) == 0 || compound(1,O_index) == 0
            is_simple = false;
        else
            if compound(1,N_index) ~= 0
                if compound(1,P_index) ~= 0 || compound(1,S_index) ~= 0
                    is_simple = false;
                else
                    % CHON
                    is_simple = true;
                end
            elseif compound(1,P_index) ~= 0
                if compound(1,S_index) ~= 0
                    is_simple = false;
                else
                    % CHOP
                    is_simple = true;
                end
            else
              % CHO & CHOS
              is_simple = true;
            end
        end
    end
end

% Formatting function for data export into Excel
function data_reformatted = FormatDataForExport(peaklist)
global format

    mz=peaklist(:,format.Column_mz);
    Magnitude=peaklist(:,format.Column_Magnitude);
    SN=peaklist(:,format.Column_SN);
    EstimC=peaklist(:,format.Column_EstimC);
    C=peaklist(:,format.Column_C);
    H=peaklist(:,format.Column_H);
    O=peaklist(:,format.Column_O);
    N=peaklist(:,format.Column_N);
    S=peaklist(:,format.Column_S);
    P=peaklist(:,format.Column_P);
    E=peaklist(:,format.Column_E);
    ExactMass=peaklist(:,format.Column_ExactMass);
    e=peaklist(:,format.Column_error);
    OC=peaklist(:,format.Column_OC);
    HC=peaklist(:,format.Column_HC);
    DBE=peaklist(:,format.Column_DBE);
    DBEC=peaklist(:,format.Column_DBEC);
    AImod=peaklist(:,format.Column_AImod);
    
    [r,~]=size(peaklist);    % Determine type of formula
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

% Reformat
data_reformatted=[num2str(mz,'%4.10f'), num2str(Magnitude), num2str(SN), num2str(EstimC), string(t), ...
    num2str(C), num2str(H), num2str(O), num2str(N), num2str(S), num2str(P), num2str(E), ... 
    num2str(ExactMass,'%4.10f'), num2str(e),num2str(OC), num2str(HC), num2str(DBE), num2str(DBEC), num2str(AImod)];
end