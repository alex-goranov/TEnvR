function FTMS_FormulaAssignment(sample,config)
%% Description

% This code generates molecular formulas for a list of m/z values. 

% Example:
% FTMS_FormulaAssignment('Sample 1_Refinement.txt',FTMS_ConfigurationAssignment)

% The algorhitm of this code has been previously published: 

% Obeid, W.A. (2015) Investigation of the potential for algaenan to produce 
% hydrocarbon based fuels from algae by hydrous pyrolysis, Department of 
% Chemistry and Biochemistry. Old Dominion University, pp. 1-223.

%% Code

tic

Data=load(sample);
Data=sortrows(Data,1);

mz = Data(:,1);             %m/z peak
magnitude = Data(:,2);      %Intensity

% Establish ranges for each element (except C)
Range_H = config.H_min:config.H_max;    
Range_O = config.O_min:config.O_max; 
Range_N = config.N_min:config.N_max; 
Range_S = config.S_min:config.S_max; 
Range_P = config.P_min:config.P_max; 
Range_E = config.E_min:config.E_max; 
Range_K = config.K_min:config.K_max; 
Range_Na = config.Na_min:config.Na_max; 

NumberOfPeaks = size(mz,1);
AllPeaks = cell(NumberOfPeaks,1);

% Creates a matrix of all possible formula combinations
if strcmpi(config.Ionization_Mode,'negative')
    FormulaCombinations = combvec(config.C_min,Range_H,Range_N,Range_O,Range_S,Range_P,Range_E); 
elseif strcmpi(config.Ionization_Mode,'positive')
    FormulaCombinations = combvec(config.C_min,Range_H,Range_N,Range_O,Range_S,Range_P,Range_E,Range_K,Range_Na);
else
    error('Error! Need to select the correct ionization mode (positive or negative) in the configuration file!')
end

% Establish a format of the elements in the formulas
format.Row_C=1;
format.Row_H=2;
format.Row_N=3;
format.Row_O=4;
format.Row_S=5;
format.Row_P=6;
format.Row_E=7;
format.Row_K=8;
format.Row_Na=9;

% Refine possible formula combinations by DBE
if strcmpi(config.Ionization_Mode,'negative')
    if config.Heteroelement_Halogen && ~config.Heteroelement_N15
         DBE = 1+(((2*FormulaCombinations(format.Row_C,:))-(FormulaCombinations(format.Row_H,:)+1)+FormulaCombinations(format.Row_N,:)+FormulaCombinations(format.Row_P,:)-FormulaCombinations(format.Row_E,:))/2);
    elseif ~config.Heteroelement_Halogen && config.Heteroelement_N15    
        DBE = 1+(((2*FormulaCombinations(format.Row_C,:))-(FormulaCombinations(format.Row_H,:)+1)+FormulaCombinations(format.Row_N,:)+FormulaCombinations(format.Row_P,:)+FormulaCombinations(format.Row_E,:))/2);
    else
        DBE = 1+(((2*FormulaCombinations(format.Row_C,:))-(FormulaCombinations(format.Row_H,:)+1)+FormulaCombinations(format.Row_N,:)+FormulaCombinations(format.Row_P,:))/2);
    end
    
    Index_Refined=find(floor(DBE) == DBE & DBE>=config.DBE_min & DBE<=config.DBE_max); % Find combinations with a whole AND positive DBE value
    FormulaCombinations_Refined = FormulaCombinations(:,Index_Refined);
    
    ExactMass = FormulaCombinations_Refined(format.Row_C,:)*config.Mass_12C + FormulaCombinations_Refined(format.Row_H,:)*config.Mass_1H + FormulaCombinations_Refined(format.Row_N,:)*config.Mass_14N +...
        FormulaCombinations_Refined(format.Row_O,:)*config.Mass_16O + FormulaCombinations_Refined(format.Row_S,:)*config.Mass_32S + FormulaCombinations_Refined(format.Row_P,:)*config.Mass_31P + ...
        FormulaCombinations_Refined(format.Row_E,:)*config.Mass_Heteroelement + config.Mass_electron; 

elseif strcmpi(config.Ionization_Mode,'positive')
    % Refine Formula Combinations by removing formulas with both Na and K
    NaK=FormulaCombinations(format.Row_K,:)+FormulaCombinations(format.Row_Na,:);
    index_refined1=find(NaK <= 1); % Find combinations with Na+K<=1
    FormulaCombinations_Refined1 = FormulaCombinations(:,index_refined1);
    
    % Refine Formula Combinations by DBE
    if config.Heteroelement_Halogen && ~config.Heteroelement_N15
        DBE = 1+(((2*FormulaCombinations_Refined1(format.Row_C,:))-(FormulaCombinations_Refined1(format.Row_H,:)-1+FormulaCombinations_Refined1(format.Row_K,:)+FormulaCombinations_Refined1(format.Row_Na,:))+FormulaCombinations_Refined1(format.Row_N,:)+FormulaCombinations_Refined1(format.Row_P,:)-FormulaCombinations_Refined1(format.Row_E,:))/2);
    elseif ~config.Heteroelement_Halogen && config.Heteroelement_N15
        DBE = 1+(((2*FormulaCombinations_Refined1(format.Row_C,:))-(FormulaCombinations_Refined1(format.Row_H,:)-1+FormulaCombinations_Refined1(format.Row_K,:)+FormulaCombinations_Refined1(format.Row_Na,:))+FormulaCombinations_Refined1(format.Row_N,:)+FormulaCombinations_Refined1(format.Row_P,:)+FormulaCombinations_Refined1(format.Row_E,:))/2);
    else
        DBE = 1+(((2*FormulaCombinations_Refined1(format.Row_C,:))-(FormulaCombinations_Refined1(format.Row_H,:)-1+FormulaCombinations_Refined1(format.Row_K,:)+FormulaCombinations_Refined1(format.Row_Na,:))+FormulaCombinations_Refined1(format.Row_N,:)+FormulaCombinations_Refined1(format.Row_P,:))/2);
    end
    
    index_refined2=find(floor(DBE) == DBE & DBE>=config.DBE_min & DBE<=config.DBE_max); %Find combinations with a whole AND positive DBE value
    FormulaCombinations_Refined = FormulaCombinations_Refined1(:,index_refined2);
    
    ExactMass = FormulaCombinations_Refined(format.Row_C,:)*config.Mass_12C + FormulaCombinations_Refined(format.Row_H,:)*config.Mass_1H + FormulaCombinations_Refined(format.Row_N,:)*config.Mass_14N +...
        FormulaCombinations_Refined(format.Row_O,:)*config.Mass_16O + FormulaCombinations_Refined(format.Row_S,:)*config.Mass_32S + FormulaCombinations_Refined(format.Row_P,:)*config.Mass_31P + ...
        FormulaCombinations_Refined(format.Row_E,:)*config.Mass_Heteroelement + FormulaCombinations_Refined(format.Row_K,:)*config.Mass_40K + FormulaCombinations_Refined(format.Row_Na,:)*config.Mass_23Na - config.Mass_electron; 
else
    error('Error! Need to select the correct ionization mode (positive or negative) in the configuration file!')
end

for i = 1:NumberOfPeaks
    C_max = floor(mz(i)/config.Mass_12C);
    
    if C_max <= max(FormulaCombinations_Refined(format.Row_C,:)) || C_max <= config.C_min
        ;       
    else
        C_min_new = max(FormulaCombinations_Refined(format.Row_C,:));
        
        if strcmpi(config.Ionization_Mode,'negative')
            FormulaCombinations_New = combvec(C_min_new + 1:C_max,Range_H,Range_N,Range_O,Range_S,Range_P,Range_E); %Must match format!!!
            
            % Refine Formula Combinations by DBE
            if config.Heteroelement_Halogen && ~config.Heteroelement_N15
                 DBE_New = 1+(((2*FormulaCombinations_New(format.Row_C,:))-(FormulaCombinations_New(format.Row_H,:)+1)+FormulaCombinations_New(format.Row_N,:)+FormulaCombinations_New(format.Row_P,:)-FormulaCombinations_New(format.Row_E,:))/2);
            elseif ~config.Heteroelement_Halogen && config.Heteroelement_N15
                 DBE_New = 1+(((2*FormulaCombinations_New(format.Row_C,:))-(FormulaCombinations_New(format.Row_H,:)+1)+FormulaCombinations_New(format.Row_N,:)+FormulaCombinations_New(format.Row_P,:)+FormulaCombinations_New(format.Row_E,:))/2);
            else
                DBE_New = 1+(((2*FormulaCombinations_New(format.Row_C,:))-(FormulaCombinations_New(format.Row_H,:)+1)+FormulaCombinations_New(format.Row_N,:)+FormulaCombinations_New(format.Row_P,:))/2);
            end
            
            Index_Refined_New=find(floor(DBE_New) == DBE_New & DBE_New>=config.DBE_min & DBE_New<=config.DBE_max); % Find combinations with a whole AND positive DBE value
            FormulaCombinations_New_Refined = FormulaCombinations_New(:,Index_Refined_New);
            
            ExactMass_New = FormulaCombinations_New_Refined(format.Row_C,:)*config.Mass_12C + FormulaCombinations_New_Refined(format.Row_H,:)*config.Mass_1H + FormulaCombinations_New_Refined(format.Row_N,:)*config.Mass_14N +...
                FormulaCombinations_New_Refined(format.Row_O,:)*config.Mass_16O + FormulaCombinations_New_Refined(format.Row_S,:)*config.Mass_32S + FormulaCombinations_New_Refined(format.Row_P,:)*config.Mass_31P + ...
                FormulaCombinations_New_Refined(format.Row_E,:)*config.Mass_Heteroelement + config.Mass_electron;      
        elseif strcmpi(config.Ionization_Mode,'positive')
            FormulaCombinations_New = combvec(C_min_new + 1:C_max,Range_H,Range_N,Range_O,Range_S,Range_P,Range_E,Range_K,Range_Na); %Must match format!!!

            % Refine Formula Combinations by removing formulas with both Na and K
            NaK_New=FormulaCombinations_New(format.Row_K,:)+FormulaCombinations_New(format.Row_Na,:);
            index_refined1_New=find(NaK_New <= 1); % Find combinations with Na+K<=1
            FormulaCombinations_New_Refined1 = FormulaCombinations_New(:,index_refined1_New);

            % Refine Formula Combinations by DBE
            if config.Heteroelement_Halogen && ~config.Heteroelement_N15
                DBE_New = 1+(((2*FormulaCombinations_New_Refined1(format.Row_C,:))-(FormulaCombinations_New_Refined1(format.Row_H,:)-1+FormulaCombinations_New_Refined1(format.Row_K,:)+FormulaCombinations_New_Refined1(format.Row_Na,:))+FormulaCombinations_New_Refined1(format.Row_N,:)+FormulaCombinations_New_Refined1(format.Row_P,:)-FormulaCombinations_New_Refined1(format.Row_E,:))/2);
            elseif ~config.Heteroelement_Halogen && config.Heteroelement_N15
                DBE_New = 1+(((2*FormulaCombinations_New_Refined1(format.Row_C,:))-(FormulaCombinations_New_Refined1(format.Row_H,:)-1+FormulaCombinations_New_Refined1(format.Row_K,:)+FormulaCombinations_New_Refined1(format.Row_Na,:))+FormulaCombinations_New_Refined1(format.Row_N,:)+FormulaCombinations_New_Refined1(format.Row_P,:)+FormulaCombinations_New_Refined1(format.Row_E,:))/2);
            else
                DBE_New = 1+(((2*FormulaCombinations_New_Refined1(format.Row_C,:))-(FormulaCombinations_New_Refined1(format.Row_H,:)-1+FormulaCombinations_New_Refined1(format.Row_K,:)+FormulaCombinations_New_Refined1(format.Row_Na,:))+FormulaCombinations_New_Refined1(format.Row_N,:)+FormulaCombinations_New_Refined1(format.Row_P,:))/2);
            end

            index_refined2_New=find(floor(DBE_New) == DBE_New & DBE_New>=config.DBE_min & DBE_New<=config.DBE_max); %Find combinations with a whole AND positive DBE value
            FormulaCombinations_New_Refined = FormulaCombinations_New_Refined1(:,index_refined2_New);

            ExactMass_New = FormulaCombinations_New_Refined(format.Row_C,:)*config.Mass_12C + FormulaCombinations_New_Refined(format.Row_H,:)*config.Mass_1H + FormulaCombinations_New_Refined(format.Row_N,:)*config.Mass_14N +...
                FormulaCombinations_New_Refined(format.Row_O,:)*config.Mass_16O + FormulaCombinations_New_Refined(format.Row_S,:)*config.Mass_32S + FormulaCombinations_New_Refined(format.Row_P,:)*config.Mass_31P + ...
                FormulaCombinations_New_Refined(format.Row_E,:)*config.Mass_Heteroelement + FormulaCombinations_New_Refined(format.Row_K,:)*config.Mass_40K + FormulaCombinations_New_Refined(format.Row_Na,:)*config.Mass_23Na - config.Mass_electron; 
        else
            error('Error! Need to select the correct ionization mode (positive or negative) in the configuration file!')
        end
        
        FormulaCombinations_Refined = [FormulaCombinations_Refined FormulaCombinations_New_Refined];
        ExactMass = [ExactMass ExactMass_New];
    end
    
    AssignmentError =((mz(i)-ExactMass)./ExactMass)*10^6;
    
    % Identify formulas within the specified Assignment Accuracy:
    location = find(AssignmentError >= -config.MassAccuracy & AssignmentError <= config.MassAccuracy);
    
    A = repmat(size(location,2),1,size(location,2))';   % A is number of possibilities
    B = repmat(mz(i),1,size(location,2))';              % B is m/z from instrument
    C = repmat(magnitude(i),1,size(location,2))';       % C is Peak Height 
    D = FormulaCombinations_Refined(:,location)';       % D identifies correct formulas
    E = ExactMass(location)';                           % E is exact mass of each formula
    F = AssignmentError(location)';                     % F is error on each formula
    
    AllPeaks{i,1} = [ A B C D E F ];  
end

Formulas_Assigned = cell2mat(AllPeaks);

%Add two empty columns for K and Na in negative-ESI. They will be already there in positive-ESI.
if strcmpi(config.Ionization_Mode,'negative')
   Formulas_Assigned=[Formulas_Assigned(:,1:10),zeros(size(Formulas_Assigned,1),2),Formulas_Assigned(:,11:12)]; 
end

T=table(Formulas_Assigned);
T = splitvars(T, 'Formulas_Assigned');
T.Properties.VariableNames{1} = 'Options';
T.Properties.VariableNames{2} = 'm/z';
T.Properties.VariableNames{3} = 'Magnitude';
T.Properties.VariableNames{4} = 'C';
T.Properties.VariableNames{5} = 'H';
T.Properties.VariableNames{6} = 'N';
T.Properties.VariableNames{7} = 'O';
T.Properties.VariableNames{8} = 'S';
T.Properties.VariableNames{9} = 'P';
T.Properties.VariableNames{10} = 'E';
T.Properties.VariableNames{11} = 'K';
T.Properties.VariableNames{12} = 'Na';
T.Properties.VariableNames{13} = 'm/z2';
T.Properties.VariableNames{14} = 'error';
         
writetable(T,[sample(1:end-4) '_F.txt'],'Delimiter','\t','WriteRowNames',true);

disp(['Finished assigning formulas to the peak list ' char([sample(1:end-4) '.txt' ' (' num2str(toc) ' seconds)'])])
end