function varargout = FTMS_RefinementPeaks(filename,blank,config)
%% Description

%This code takes an FT-ICR-MS file with calibrated data for a sample and a blank and:
    %1) Identifies blank peaks (found in the blank file)
    %2) Identifies inorganic "salt" peaks (see below)
    %3) Identified isotopologue peaks (see below)
    %4) Prepares a peak list of peaks ready for formula assignment

% Blank peaks: peaks present in both sample and blank (within 0.0005 m/z)
% Inorganic (salt) peaks: mass with a mass defect of %0.4-0.99 up to 400m/z and 0.6-0.99 for 400-800)
% 13C isotopologue peaks: 13C:12C peak ratio < 0.5; C number < 46 (i.e. 0.5/0.011);
% 34S isotopologue peaks: 34S:32S peak ratio < 0.3; S number < 6.6 (i.e. 0.3/0.045);
% 54Fe isotopologue peaks: 54Fe:56Fe peak ratio < 0.2 Fe number < 3.2 (i.e. 0.2/0.063); 
% 37Cl isotopologue peaks: 37Cl:35Cl peak ratio 0.2 - 2; Cl number 0.6 - 6.2 (i.e. 2/0.32);
% 200Hg isotopologue peaks: 202Hg peak ratio 0.5 - 2; Hg number 0.65 - 2.6 (i.e. 2/0.77);

% The sample and blank peak lists must be located in text (.txt) files with
% no text headers!

% Portions of this code are obtained from previously published codes: 

% Obeid, W.A. (2015) Investigation of the Potential for Algaenan to Produce Hydrocarbon 
% Based Fuels From Algae by Hydrous Pyrolysis, Department of Chemistry and Biochemistry. 
% Old Dominion University, pp. 1-223.

% Patriarca, C. and Hawkes, J.A. (2020) High Molecular Weight Spectral Interferences in 
% Mass Spectra of Dissolved Organic Matter. Journal of the American Society for Mass Spectrometry 32, 394-397.

% Examples:
% FTMS_RefinementPeaks('Sample 1.txt','Blank PPL.txt',FTMS_ConfigurationAssignment)

%% Stage 1: Data Import and trimming

tic
dataSample_raw = load(filename);  % load mass list of your sample from a text data file
dataBlank_raw = load(blank);      % load mass list of the blank from a text data file

% Create ghost S/N values if they are not loaded as third column
size_dataSample_raw=size(dataSample_raw);
size_dataBlank_raw=size(dataBlank_raw);

if size_dataSample_raw(1,2) == 2
    dataSample_raw=[dataSample_raw,zeros(size_dataSample_raw(1,1),1)];
elseif size_dataSample_raw(1,2) == 3
    ;
else
    disp('Warning! Columns 4 and beyond of the loaded mass list will not be imported!')
end

if size_dataBlank_raw(1,2) == 2
    dataBlank_raw=[dataBlank_raw,zeros(size_dataBlank_raw(1,1),1)];
elseif size_dataSample_raw(1,2) == 3
    ;
else
    disp('Warning! Columns 4 and beyond of the loaded mass list will not be imported!')
end

dataSample_raw=sortrows(dataSample_raw,1);
dataBlank_raw=sortrows(dataBlank_raw,1);

% Trimming the mass spectra to based on previously defined cutoffs
indeces_sample = dataSample_raw(:,1)>=config.MZcutoff_low & dataSample_raw(:,1)<=config.MZcutoff_high;
dataSample=dataSample_raw(indeces_sample,:);

indeces_blank = dataBlank_raw(:,1)>=config.MZcutoff_low+1 & dataBlank_raw(:,1)<=config.MZcutoff_high+1;
dataBlank=dataBlank_raw(indeces_blank,:);

% Extracting the data
SampleMZ=dataSample(:,1);    % m/z values of the SAMPLE, column 1
SampleInt=dataSample(:,2);   % peak magnitudes of the SAMPLE, column 2
SampleSN=dataSample(:,3);    % S/N values of the SAMPLE, column 3

BlankMZ=dataBlank(:,1);      % m/z values of the BLANK, column 1
BlankInt=dataBlank(:,2);     % peak magnitudes of the BLANK, column 2
BlankSN=dataBlank(:,3);      % S/N values of the BLANK, column 3

%% Stage 2: Identify blank peaks:

%C ompare m/z of sample and blank peaks:
PeaksBlank = zeros(size(SampleMZ,1), 2);
for i = 1 : size(SampleMZ,1)
    for j = 2 : size(BlankMZ,1)
        p = SampleMZ(i) - BlankMZ(j);

        if abs(((SampleMZ(i)-BlankMZ(j))/SampleMZ(i))*1e6) < config.MassAccuracy
            PeaksBlank(i, 1) = BlankMZ(j);
            PeaksBlank(i, 2) = BlankInt(j);
            PeaksBlank(i, 3) = BlankSN(j);
        end
    end
end

%% Stage 3: Identify inorganic (salt) peaks:

PeaksSalt = zeros(size(SampleMZ,1), 2);
for i = 1 : size(SampleMZ,1)
    Massdef(i)=SampleMZ(i)-fix(SampleMZ(i)); % Calculate mass defect
    if SampleMZ(i) < 300 && Massdef(i) >= 0.4 && Massdef(i)<=0.97
       PeaksSalt(i, 1) = SampleMZ(i);
       PeaksSalt(i, 2) = SampleInt(i);
       PeaksSalt(i, 3) = SampleSN(i);
    elseif SampleMZ(i)>= 300 && SampleMZ(i)<=500 && Massdef(i)>=0.5 && Massdef(i)<=0.96
       PeaksSalt(i, 1) = SampleMZ(i);
       PeaksSalt(i, 2) = SampleInt(i);
       PeaksSalt(i, 3) = SampleSN(i);
   elseif SampleMZ(i)>= 500 && SampleMZ(i)<=800 && Massdef(i)>=0.6 && Massdef(i)<=0.95
       PeaksSalt(i, 1) = SampleMZ(i);
       PeaksSalt(i, 2) = SampleInt(i);
       PeaksSalt(i, 3) = SampleSN(i);
    elseif SampleMZ(i) > 800  && Massdef(i)>=0.7 && Massdef(i)<=0.94
       PeaksSalt(i, 1) = SampleMZ(i);
       PeaksSalt(i, 2) = SampleInt(i);
       PeaksSalt(i, 3) = SampleSN(i);
    end
end

%% Stage 4: Identify 13C isotopologues

% identify all the 12C & 13C by mass: 1.003355
iso13C = zeros(size(SampleMZ,1), 2);
for i = 1 :(size(SampleMZ,1)-1)
    for j = 2 : size(SampleMZ,1)
        if j < i
        continue;
        else 
    p = SampleMZ(j) - SampleMZ(i);
    a = SampleInt(j) - 0.5 * SampleInt(i); % 50% threshold 
    
    if abs(((SampleMZ(j)-1.003355-SampleMZ(i))/SampleMZ(i))*1e6) < config.MassAccuracy && a < 0
        iso13C(i, 2) = (SampleInt(j)/SampleInt(i))*100/1.1; 
        iso13C(i, 3) = SampleMZ(j);
        iso13C(i, 4) = SampleInt(j);
        iso13C(i, 5) = p;
        iso13C(j, 1) = SampleMZ(j);       
    end
        end
    end
   
end

if sum(iso13C(:,1)) == 0
    iso13C=[iso13C,zeros(size(iso13C(:,1),1),3)];
end

%% Stage 5: Identify 34S isotopologues

% identify all the 32S & 34S by mass: 1.995796
iso34S = zeros(size(SampleMZ,1), 2);
for i = 1 :(size(SampleMZ,1)-1)
    for j = 2 : size(SampleMZ,1)
        if j < i
        continue;
        else 
    p = SampleMZ(j) - SampleMZ(i);
    a = SampleInt(j) - 0.3 * SampleInt(i); % 30% threshold
    
    if abs(((SampleMZ(j)-1.995796-SampleMZ(i))/SampleMZ(i))*1e6) < config.MassAccuracy && a < 0
        iso34S(i, 2) = (SampleInt(j)/SampleInt(i))*100/4.5;
        iso34S(i, 3) = SampleMZ(j);
        iso34S(i, 4) = SampleInt(j);
        iso34S(i, 5) = p;
        iso34S(j, 1) = SampleMZ(j);       
    end
        end
    end
   
end

if sum(iso34S(:,1)) == 0
    iso34S=[iso34S,zeros(size(iso34S(:,1),1),3)];
end

%% Stage 6: Identify 54Fe isotopologues

% identify all the 56Fe & 54Fe by mass: 1.995327
iso54Fe = zeros(size(SampleMZ,1), 2);
for i = 2 :size(SampleMZ,1)
    for j = 1 : (size(SampleMZ,1)-1)
        if j > i
        continue;
        else 
    p = SampleMZ(i) - SampleMZ(j);
    a = SampleInt(j) - 0.2 * SampleInt(i); % 20% threshold
    
    if abs(((SampleMZ(i)-1.995327-SampleMZ(j))/SampleMZ(j))*1e6) < config.MassAccuracy && a < 0
        iso54Fe(i, 2) = (SampleInt(j)/SampleInt(i))*100/6.3;
        iso54Fe(i, 3) = SampleMZ(j);
        iso54Fe(i, 4) = SampleInt(j);
        iso54Fe(i, 5) = p;
        iso54Fe(j, 1) = SampleMZ(j);       
    end
        end
    end
   
end

if sum(iso54Fe(:,1)) == 0
    iso54Fe=[iso54Fe,zeros(size(iso54Fe(:,1),1),3)];
end

%% Stage 7: Identify 37Cl isotopologues

% identify all the 35Cl 37Cl by mass: 1.99705
iso37Cl = zeros(size(SampleMZ,1), 2);
for i = 1 :(size(SampleMZ,1)-1)
    for j = 2 : size(SampleMZ,1)
        if j < i
        continue;
        else 
    p = SampleMZ(j) - SampleMZ(i);
    a = SampleInt(j) - 0.2 * SampleInt(i);
    b = SampleInt(j) - 2 * SampleInt(i);
    
    if abs(((SampleMZ(j)-1.99705-SampleMZ(i))/SampleMZ(i))*1e6) < config.MassAccuracy && a > 0 && b<0
        iso37Cl(i, 2) = (SampleInt(j)/SampleInt(i))*100/32;
        iso37Cl(i, 3) = SampleMZ(j);
        iso37Cl(i, 4) = SampleInt(j);
        iso37Cl(i, 5) = p;
        iso37Cl(j, 1) = SampleMZ(j);        
    end
        end
    end
   
end

if sum(iso37Cl(:,1)) == 0
    iso37Cl=[iso37Cl,zeros(size(iso37Cl(:,1),1),3)];
end

%% Stage 8: Identify 202Hg isotopologues

% identify all the 202Hg 200Hg by mass: 2.002316
iso200Hg = zeros(size(SampleMZ,1), 2);
for i = 2 :size(SampleMZ,1)
    for j = 1 : (size(SampleMZ,1)-1)
        if j > i
        continue;
        else 
     p = SampleMZ(i) - SampleMZ(j);
     a = SampleInt(j) - 0.5 * SampleInt(i);
     b = SampleInt(j) - 2 * SampleInt(i);
     
    if abs(((SampleMZ(i)-2.002316-SampleMZ(j))/SampleMZ(j))*1e6) < config.MassAccuracy && a > 0 && b < 0
        iso200Hg(i, 2) = (SampleInt(j)/SampleInt(i))*100/77;
        iso200Hg(i, 3) = SampleMZ(j);
        iso200Hg(i, 4) = SampleInt(j);
        iso200Hg(i, 5) = p;
        iso200Hg(j, 1) = SampleMZ(j);       
    end
        end
    end
   
end

if sum(iso200Hg(:,1)) == 0
    iso200Hg=[iso200Hg,zeros(size(iso200Hg(:,1),1),3)];
end

%% Stage 10: Identify doubly-charged ions

PeaksDoublyCharged=zeros(size(SampleMZ,1),2);
for m=1:size(SampleMZ,1)                        % SampleMZ being the 13C isotopologue
    monoisotopic_mass=SampleMZ(m)-1.003355/2;   % Finging out the 12C monoisotopic peak
    m_diff=abs(SampleMZ-monoisotopic_mass);   
    monoisotopic_position=find(m_diff==min(m_diff)); % Finding where is the position of the 12C monoisotopic peak is
    a = SampleInt(m) - 0.5 * SampleInt(monoisotopic_position);
    for k=1:size(monoisotopic_position,1)
        monoisotopic_position_k=monoisotopic_position(k);
        if abs(((SampleMZ(monoisotopic_position_k)-monoisotopic_mass)/monoisotopic_mass)*1e6) < config.MassAccuracy && a < 0 && PeaksSalt(monoisotopic_position_k) == 0        
            PeaksDoublyCharged(monoisotopic_position_k,1)=SampleMZ(monoisotopic_position_k);
            PeaksDoublyCharged(monoisotopic_position_k,2)=SampleInt(monoisotopic_position_k);
            PeaksDoublyCharged(m,1)=SampleMZ(m);
            PeaksDoublyCharged(m,2)=SampleInt(m);
            PeaksDoublyCharged(m,3)=SampleSN(m);
        end  
    end
end

%% Stage 9: Export an Excel report with data and statistics,as well as a "refined" peak list

warning('off','MATLAB:xlswrite:AddSheet')
DataAll = [dataSample, PeaksBlank(:,1), PeaksSalt(:,1), PeaksDoublyCharged(:,1), iso13C, iso34S, iso54Fe, iso37Cl, iso200Hg];
Titles={'m/z', 'Magnitude', 'S/N', 'Blank', 'Salt','Doubly Charged', '13C','Estim C#','13C m/z','13C Int','difC13','34S','Estim S#','34S m/z','34S Int','difS34','54Fe','Estim Fe#','54Fe m/z','54Fe Int','difFe','37Cl','Estim Cl#','37Cl m/z','37Cl Int','difCl','200Hg','Estim Hg#','200Hg m/z','200Hg Int','difHg'};
xlswrite([filename(1:end-4) '_Refinement.xlsx'], Titles, 1, 'A1');
xlswrite([filename(1:end-4) '_Refinement.xlsx'], DataAll, 1, 'A2');

% Remove blank, salt, doubly-charged, and isotopoluges peaks. 
% In the exported file, it will contain m/z, intensity, estimated C number, and S/N.

DataRefined = DataAll(DataAll(:,4) == 0,[1:31]);          % Removing Blank Peaks (Column 4) 
DataRefined = DataRefined(DataRefined(:,5) == 0,[1:31]);  % Removing Salt Peaks (Column 5) 
DataRefined = DataRefined(DataRefined(:,6) == 0,[1:31]);  % Removing Doubly Charged Peaks (Column 6) 
DataRefined = DataRefined(DataRefined(:,7) == 0,[1:31]);  % Removing 13C Peaks (Column 7) 
DataRefined = DataRefined(DataRefined(:,12) == 0,[1:31]); % Removing 34S Peaks (Column 12) 
DataRefined = DataRefined(DataRefined(:,17) == 0,[1:31]); % Removing 54Fe Peaks (Column 17) 
DataRefined = DataRefined(DataRefined(:,22) == 0,[1:31]); % Removing 37Cl Peaks (Column 22) 
DataRefined = DataRefined(DataRefined(:,27) == 0,[1:31]); % Removing 200Hg Peaks (Column 27) 


DataRefined_trimed=[DataRefined(:,1), DataRefined(:,2), DataRefined(:,3),DataRefined(:,8)];
DataRefined_Titles={'m/z','Magnitude','S/N','Estim C#'};
DataRefined_ExportTXT=[DataRefined(:,1),DataRefined(:,2)];

if config.Export_Cl
    DataRefinedCl = DataRefined(DataRefined(:,22) ~= 0,[1,2]);

    fileID = fopen([filename(1:end-4) '_Refinement_Cl.txt'], 'w');
    fprintf(fileID,'%.10f %.0f\r\n', DataRefinedCl');
    fclose(fileID);   
end

fileID = fopen([filename(1:end-4) '_Refinement.txt'], 'w');
fprintf(fileID,'%.10f %.0f\r\n', DataRefined_ExportTXT');
fclose(fileID);

xlswrite([filename(1:end-4) '_Refinement.xlsx'], DataRefined_Titles, 'DataRefined', 'A1');
xlswrite([filename(1:end-4) '_Refinement.xlsx'], DataRefined_trimed, 'DataRefined', 'A2');

%% Stage 9: Statistics

[Stats_Num_AllPeaks,~]=size(SampleInt);
[Stats_Num_BlankPeaks,~]=size(PeaksBlank(PeaksBlank(:,2) ~=0));
[Stats_Num_SaltPeaks,~]=size(PeaksSalt(PeaksSalt(:,2) ~=0));
[Stats_Num_DoublyChargedPeaks,~]=size(PeaksDoublyCharged(PeaksDoublyCharged(:,2) ~=0));
[Stats_Num_13CPeaks,~]=size(iso13C(iso13C(:,1) ~=0));
[Stats_Num_34SPeaks,~]=size(iso34S(iso34S(:,1) ~=0));
[Stats_Num_54FePeaks,~]=size(iso54Fe(iso54Fe(:,1) ~=0));
[Stats_Num_37ClPeaks,~]=size(iso37Cl(iso37Cl(:,1) ~=0));
[Stats_Num_200HgPeaks,~]=size(iso200Hg(iso200Hg(:,1) ~=0));
[Stats_Num_RefinedPeaks,~]=size(DataRefined(DataRefined(:,1) ~=0));
Stats_Num_BadPeaks=Stats_Num_AllPeaks-Stats_Num_RefinedPeaks;

Stats_Int_AllPeaks=sum(SampleInt);
Stats_Int_BlankPeaks=sum(PeaksBlank(:,2));
Stats_Int_SaltPeaks=sum(PeaksSalt(:,2));
Stats_Int_DoublyChargedPeaks=sum(PeaksDoublyCharged(:,2));
Stats_Int_13CPeaks=sum(iso13C(:,4));
Stats_Int_34SPeaks=sum(iso34S(:,4));
Stats_Int_54FePeaks=sum(iso54Fe(:,4));
Stats_Int_37ClPeaks=sum(iso37Cl(:,4));
Stats_Int_200HgPeaks=sum(iso200Hg(:,4));
Stats_Int_RefinedPeaks=sum(DataRefined(:,2));
Stats_Int_BadPeaks=Stats_Int_AllPeaks-Stats_Int_RefinedPeaks;

% Statistics Matrix
Stats_TitlesRows={'All Peaks';'Blank'; 'Salt';'Doubly Charged';'13C';'34S';'54Fe';'37Cl';'200Hg';'Rejected Peaks';'Refined Peaks'};
Stats_TitlesColumns={'', '# of peaks', '% Num', '% Magn'};
Stats_note={'Please note that # of bad peaks may be < Blank+Salt+isotope peaks due to overlaping peaks (i.e. a blank peak may also be a salt peak!)'};
Stats=[Stats_Num_AllPeaks,100,100;...
    Stats_Num_BlankPeaks,round((Stats_Num_BlankPeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_BlankPeaks*100/Stats_Int_AllPeaks),2);...
    Stats_Num_SaltPeaks,round((Stats_Num_SaltPeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_SaltPeaks*100/Stats_Int_AllPeaks),2);...
    Stats_Num_DoublyChargedPeaks,round((Stats_Num_DoublyChargedPeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_DoublyChargedPeaks*100/Stats_Int_AllPeaks),2);...
    Stats_Num_13CPeaks,round((Stats_Num_13CPeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_13CPeaks*100/Stats_Int_AllPeaks),2);...
    Stats_Num_34SPeaks,round((Stats_Num_34SPeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_34SPeaks*100/Stats_Int_AllPeaks),2);...
    Stats_Num_54FePeaks,round((Stats_Num_54FePeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_54FePeaks*100/Stats_Int_AllPeaks),2);...
    Stats_Num_37ClPeaks,round((Stats_Num_37ClPeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_37ClPeaks*100/Stats_Int_AllPeaks),2);...
    Stats_Num_200HgPeaks,round((Stats_Num_200HgPeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_200HgPeaks*100/Stats_Int_AllPeaks),2);...
    Stats_Num_BadPeaks,round((Stats_Num_BadPeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_BadPeaks*100/Stats_Int_AllPeaks),2);...
    Stats_Num_RefinedPeaks,round((Stats_Num_RefinedPeaks*100/Stats_Num_AllPeaks),2),round((Stats_Int_RefinedPeaks*100/Stats_Int_AllPeaks),2)];
Stats=string(Stats);  

xlswrite([filename(1:end-4) '_Refinement.xlsx'], Stats_TitlesColumns, 1, 'AG2');
xlswrite([filename(1:end-4) '_Refinement.xlsx'], Stats_TitlesRows, 1, 'AG3');
xlswrite([filename(1:end-4) '_Refinement.xlsx'], Stats, 1, 'AH3');
xlswrite([filename(1:end-4) '_Refinement.xlsx'], Stats_note, 1, 'AL12');

if nargout == 1
    varargout{1}=[filename(1:end-4);Stats(:,1)];
else
    ;
end
disp(['Finished refining the peak list ' char(filename) ' (' num2str(toc) ' seconds)'])
end