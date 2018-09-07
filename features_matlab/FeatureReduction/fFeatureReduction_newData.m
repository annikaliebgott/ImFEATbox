function TransformedData = fFeatureReduction_newData(NormSamples, Input)
% This script performs the mapping of new data to a reduced feature space
% which has been previously computed with similar data.
% NOTE: there is no explicit mapping for every reduction algorithm. For the
% reduction techniques where there is no explicit mapping of new data to
% the reduced feature space available, the mapping is estimated from the
% orginal feature matrix and the calculated reduced feature matrix.
%
% Input variables:
%   NormSamples: normalized dataset (feature matrix)
%   Input: struct with fields:
%       sAlgo: which feature reduction algorithm should be used
%       mapping: mapping information of the previously computated
%                transformation (if available)
%       OriginalData: dataset on which the feature reduction has been
%                     performed originally
%       OriginalMappedX: orignially calculated reduced feature matrix
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic
% and Interventional Radiology, University Hospital of Tuebingen, Germany
% and the Institute of Signal Processing and System Theory University of
% Stuttgart, Germany. Last modified: September 2018
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
% ************************************************************************
%
% Implementation based on:  
%       (C) Laurens van der Maaten, Delft University of Technology
%       (c) Deng Cai, Zhejiang University

sAlgo = Input.sAlgo;

switch sAlgo   
    case {'PCA','LDA'}
        mapping = Input.mapping;
        TransformedData = bsxfun(@minus, NormSamples, mapping.mean) * mapping.M;
    case {'DiffusionMaps','Sammon','MDS','SNE','SPE','GPLVM'}
        try 
           OriginalData = Input.OriginalData;
        catch
           disp('Error: there is no explicit calculation for the reduction of new data with the algorithm(s) you specified.')
           disp('To estimate the reduction, you need to provide the original dataset and the orignal calculated mapping.')
           return
        end    
         try 
           OriginalMappedX = Input.OriginalMappedX;
        catch
           disp('Error: there is no explicit calculation for the reduction of new data with the algorithm(s) you specified.')
           disp('To estimate the reduction, you need to provide the original dataset and the original calculated mapping.')
           return
        end    
        TransformedData = fEstimateReduction_newData(NormSamples, OriginalData, OriginalMappedX);
    case 'Isomap'
        mapping = Input.mapping;
        % Precomputations for speed
        invVal = inv(diag(mapping.val));
        [val, index] = sort(mapping.val, 'descend');
        mapping.landmarks = 1:size(mapping.X, 1);
        
        val = val(1:mapping.no_dims);
        meanD1 = mean(mapping.DD .^ 2, 1);
        meanD2 = mean(mean(mapping.DD .^ 2));
        
        TransformedData = zeros(size(NormSamples, 1), mapping.no_dims);
        for i=1:size(NormSamples, 1)
            
            % Compute distance of new sample to training points
            point = NormSamples(i,:);
            tD = L2_distance(point', mapping.X');
            [~, ind] = sort(tD);
            tD(ind(mapping.k + 2:end)) = 0;
            tD = sparse(tD);
            tD = dijkstra([0 tD; tD' mapping.D], 1);
            tD = tD(mapping.landmarks + 1) .^ 2;
            
            % Compute point embedding
            subB = -.5 * (bsxfun(@minus, tD, mean(tD, 2)) - meanD1 - meanD2);
            
            vec = subB * mapping.vec * invVal;
            vec = vec(:,index(1:mapping.no_dims));
            
            TransformedData(i,:) = real(vec .* sqrt(val)');
        end
        
    case {'IsoP','OLPP'}    %matlabfunc isop algo
        eigvector = Input.eigvector;
        TransformedData = NormSamples*eigvector;        
end

end
