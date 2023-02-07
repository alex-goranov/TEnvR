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

classdef FTMS_CandidateCollection < handle
    properties
        CandidateMap
    end

    methods
        function obj = FTMS_CandidateCollection(compound_list,config)
            global format

            if nargin ~= 0
                elementMaps = num2cell(compound_list(:, format.Column_C:format.Column_C+size(format.elements,2)-1) ~= 0,2);
                composition = cellfun(@(c) char(join(format.elements(c),"",2)), elementMaps, 'UniformOutput', false);
            
                km = compound_list(:,format.Column_ExactMass)*config.Filter_KMD_Series;
                kmd_values = km - floor(km);

                compound_ids = (1:size(compound_list,1)).';
                m = size(compound_ids,1);
                counts = containers.Map("KeyType", "char", "ValueType", "double");
                for i = 1:m
                    key = composition{i};
                    if ~counts.isKey(key)
                        counts(key) = 0;
                    end
                    counts(key) = counts(key) + 1;
                end

                obj.CandidateMap = containers.Map("KeyType", "char", "ValueType", "any");
                for key = counts.keys
                    charkey = char(key);
                    candidate_mat = FTMS_Candidate();
                    candidate_mat(size(kmd_values,2), counts(charkey)) = FTMS_Candidate();
                    obj.CandidateMap(charkey) = candidate_mat;
                end

                keys = counts.keys;
                [values{1:size(keys,2)}] = deal(1);
                nextIndex = containers.Map(keys, values);
                for i = 1:m
                    key = composition{i};
                    candidate_mat = obj.CandidateMap(key);
                    for j = 1:size(kmd_values, 2)
                        candidate_mat(j,nextIndex(key)).Assign(compound_ids(i), kmd_values(i,j));
                    end
                    nextIndex(key) = nextIndex(key) + 1;
                end

                for key = obj.CandidateMap.keys
                    charkey = char(key);
                    candidate_mat = obj.CandidateMap(charkey);
                    for i = 1:size(kmd_values, 2)
                        [~, ind] = sort([candidate_mat(i,:).KMD_Value]);
                        candidate_mat(i,:) = candidate_mat(i,ind);
                    end
                    obj.CandidateMap(charkey) = candidate_mat;
                end
            end
        end

        function removed_candidates = RemoveCandidates(obj,candidates)
            candidates = sort(candidates);
            removed_candidates = containers.Map("KeyType","char","ValueType", "any");
            for key = obj.CandidateMap.keys
                charkey = char(key);
                candidate_mat = obj.CandidateMap(charkey);
                new_mat = {};
                removed_kmd_mat = {};
                for i = 1:size(candidate_mat, 1)
                    present = sort([candidate_mat(i,:).CompoundId]);
                    [~, ind] = sort([candidate_mat(i,:).CompoundId]);
                    found = ismember(present, candidates.');
                    % ind(1,~found) is the indices in original row that
                    % should be kept. sort keeps the sorting by KMD from
                    % before.
                    new_mat{i,1} = candidate_mat(i,sort(ind(1,~found)));
                    removed_kmd_mat{i,1} = candidate_mat(i,sort(ind(1,found)));
                end
                obj.CandidateMap(charkey) = vertcat(new_mat{:});
                removed_candidates(charkey) = vertcat(removed_kmd_mat{:});
            end
        end
    end
end