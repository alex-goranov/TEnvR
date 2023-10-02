classdef FTMS_KMD_Collection
    
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
    
    properties
        KMD_Map
    end

    methods
        function obj = FTMS_KMD_Collection(candidate_collection,config)
            if nargin ~= 0
                keys = candidate_collection.CandidateMap.keys;
                [values{1:size(keys,2)}] = deal(cell(size(config.Filter_KMD_Series,2),1));
                obj.KMD_Map = containers.Map(keys, values);
                for key = keys
                    charkey = char(key);
                    candidate_mat = candidate_collection.CandidateMap(charkey);
                    kmd_mat = obj.KMD_Map(charkey);
                    for i = 1:size(config.Filter_KMD_Series,2)
                        curr_kmd = FTMS_KMD_Value(candidate_mat(i,1).KMD_Value);
                        kmd_mat{i,1} = curr_kmd;
                        for j = 2:size(candidate_mat, 2)
                            candidate = candidate_mat(i, j);
                            if abs(curr_kmd.KMD - candidate.KMD_Value) < 1/10^config.Precision
                                curr_kmd.Increment()
                            else
                                curr_kmd = FTMS_KMD_Value(candidate.KMD_Value);
                                kmd_mat{i} = [kmd_mat{i}, curr_kmd];
                            end
                        end
                    end
                    obj.KMD_Map(charkey) = kmd_mat;
                end
            end
        end

        function compound_ids = FindMatchingCandidates(kmd_collection, search_values, config)
            compound_ids = [];
            for key = kmd_collection.KMD_Map.keys
                charkey = char(key);
                KMD_mat = kmd_collection.KMD_Map(charkey);
                candidate_mat = search_values.CandidateMap(charkey);
                for k=1:size(config.Filter_KMD_Series,2)
                    % i tracks kmd; j tracks search_values
                    i=1; j=1;
                    while i < size(KMD_mat{k},2) && j < size(candidate_mat, 2)
                        left = KMD_mat{k}(i).KMD;
                        right = candidate_mat(k,j).KMD_Value;
                        if abs(left-right) < 1/10^config.Precision && KMD_mat{k}(i).Count >= config.Filter_KMD_Series_Threshold
                            compound_ids = [compound_ids; candidate_mat(k,j).CompoundId];
                            j = j + 1;
                        elseif left < right
                            i = i + 1;
                        else
                            j = j + 1;
                        end
                    end
                end
            end
            compound_ids = unique(compound_ids);
        end

        function UpdateKMDValues(kmd_collection, removed_candidates, config)
            for key = kmd_collection.KMD_Map.keys
                charkey = char(key);
                KMD_mat = kmd_collection.KMD_Map(charkey);
                removed_mat = removed_candidates(charkey);
                for k=1:size(config.Filter_KMD_Series,2)
                    % i tracks collection; j tracks new values
                    i=1; j=1;
                    first = true;
                    additional_kmd = FTMS_KMD_Value.empty();
                    while i < size(KMD_mat{k},2) && j < size(removed_mat, 2)
                        left = KMD_mat{k}(i).KMD;
                        right = removed_mat(k,j).KMD_Value;
    
                        if abs(left-right) < 1/10^config.Precision
                            KMD_mat{k}(i).Increment();
                            j = j + 1;
                        elseif left < right
                            i = i + 1;
                        else
                            if first
                                to_add = FTMS_KMD_Value(right);
                                first = false;
                            elseif abs(to_add.KMD - right) >= 1/10^config.Precision
                                additional_kmd = [additional_kmd, to_add];
                                to_add = FTMS_KMD_Value(right);
                            else
                                to_add.Increment();
                            end
                            j = j + 1;
                        end
                    end
                    KMD_mat{k} = [KMD_mat{k}, additional_kmd];
                    [~, ind] = sort([KMD_mat{k}.KMD]);
                    KMD_mat{k} = KMD_mat{k}(ind);
                end
                kmd_collection.KMD_Map(charkey) = KMD_mat;
            end
        end
    end
end