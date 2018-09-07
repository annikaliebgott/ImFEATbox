function Output = fFeatureReduction(NormSample,Input)
% This Function performs the feature reduction
% Input variables:
%   NormSample: Normalized Dataset (feature matrix)
%   Input: struct with fileds:
%       Labels: Labels to Normalized Dataset
%       sAlgo: which Feature Reduction Algorithm should be used
%       N_feat_max: how many feature components should the Dataset be shrinked to
%       NN: Nearest Neighbour value. Only needed for Isomap
% Output:
%   struct with fields:
%       Transformed data: the reduced feature matrix
%       mapping: mapping information (if available)
%       eigvector: calculated eigenvectors (if available)
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


sAlgo = Input.sAlgo;
N_feat_max = Input.N_feat_max;

if strcmp(sAlgo, 'Isomap')
    try 
        NN = Input.NN;
    catch
        if ~exist('Input.NN','var')
            disp('Warning: paramter NN not defined. Set to default value NN = 10.')
            NN = 10;
        end    
    end    
end   

if any(strcmp(sAlgo,'IsoP')) || any(strcmp(sAlgo,'OLPP')) || any(strcmp(sAlgo,'LDA'))
    try
        Labels = Input.Labels;
    catch
        disp('Error: to use supervised reduction algorithms, class labels must be provided!') 
        return 
    end    
end    

switch sAlgo
    case {'PCA','DiffusionMaps','Sammon','MDS','SNE','SPE','GPLVM'}
        [mappedX, mapping] = fComputeMapping(NormSample, sAlgo, N_feat_max); 
        Output.TransformedData = mappedX;
        Output.mapping = mapping;
    case 'LDA'
        [mappedX, mapping] = lda(NormSample, Labels, N_feat_max); 
        Output.TransformedData = mappedX;
        Output.mapping = mapping;          
    case 'Isomap'
        [mappedX, mapping] = isomap(NormSample,N_feat_max,NN);
        Output.TransformedData = mappedX;
        Output.mapping = mapping;
    case 'IsoP'    
        options.k = 0;
        options.NeighborMode = 'Supervised';
        options.gnd = Labels;
        options.ReducedDim=N_feat_max;
        [eigvector, ~] = IsoP(options, NormSample);
        TransformedData = NormSample*eigvector; 
        Output.TransformedData = TransformedData;
        Output.eigvector = eigvector;
    case 'OLPP' 
        options.Metric = 'Euclidean';
        options.NeighborMode = 'Supervised';
        options.gnd = Labels;
        options.bLDA = 1;
        W = constructW(NormSample,options);
        options.PCARatio = 1;
        options.ReducedDim = N_feat_max;
        bSuccess = 0;
        while ~bSuccess
            [eigvector,~, bSuccess] = OLPP(W, options, NormSample);
        end
        TransformedData = NormSample*eigvector; 
        Output.TransformedData = TransformedData;
        Output.eigvector = eigvector;
end

end
