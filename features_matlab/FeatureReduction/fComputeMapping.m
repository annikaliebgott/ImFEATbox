function [mappedA, mapping] = fComputeMapping(A, type, no_dims, varargin)
% This function computes the mapping to the reduced feature space.
%
% ************************************************************************
% Modified for MRI feature extraction by the Department of Diagnostic
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
% Implementation by: (C) Laurens van der Maaten, 
%                        Delft University of Technology

% Check inputs
if nargin < 2
    error('Function requires at least two inputs.');
end
if ~exist('no_dims', 'var')
    no_dims = 2;
end

mapping = struct;

% Make sure there are no duplicates in the dataset
A = double(A);

% Check whether value of no_dims is correct
if ~isnumeric(no_dims) || no_dims > size(A, 2) || ((no_dims < 1 || round(no_dims) ~= no_dims) && ~any(strcmpi(type, {'PCA', 'KLM'})))
    error('Value of no_dims should be a positive integer smaller than the original data dimensionality.');
end

% Switch case
switch type
    case 'Isomap'
        % Compute Isomap mapping
        if isempty(varargin)
            [mappedA, mapping] = isomap(A, no_dims, 12);
        else
            [mappedA, mapping] = isomap(A, no_dims, varargin{1});
        end
        mapping.name = 'Isomap';
    case 'GPLVM'
        % Compute GPLVM mapping
        if isempty(varargin)
            mappedA = gplvm(A, no_dims, 1);
        else
            mappedA = gplvm(A, no_dims, varargin{1}); 
        end
        mapping.name = 'GPLVM';
    case 'DiffusionMaps'
        % Compute diffusion maps mapping
        if isempty(varargin) 
            mappedA = diffusion_maps(A, no_dims, 1, 1);
        elseif length(varargin) == 1 
            mappedA = diffusion_maps(A, no_dims, varargin{1}, 1);
        else
            mappedA = diffusion_maps(A, no_dims, varargin{1}, varargin{2}); 
        end
        mapping.name = 'DM';
    case 'SPE'
        % Compute SPE mapping
        if isempty(varargin) 
            mappedA = spe(A, no_dims, 'Global');
        elseif length(varargin) == 1
            mappedA = spe(A, no_dims, varargin{1});
        elseif length(varargin) == 2
            mappedA = spe(A, no_dims, varargin{1}, varargin{2});
        end
        mapping.name = 'SPE';
    case 'SNE'
        % Compute SNE mapping
        if isempty(varargin) 
            mappedA = sne(A, no_dims);
        else
            mappedA = sne(A, no_dims, varargin{1}); 
        end
        mapping.name = 'SNE';
    case 'LDA'
        % Run LDA on labeled dataset
        [mappedA, mapping] = lda(A(:,2:end), A(:,1), no_dims);
        mapping.name = 'LDA';
    case 'MDS'
        % Perform MDS
        mappedA = mds(A, no_dims);
        mapping.name = 'MDS';
    case 'Sammon'
        mappedA = sammon(A, no_dims);
        mapping.name = 'Sammon';
    case 'PCA'
        % Compute PCA mapping
        [mappedA, mapping] = pca(A, no_dims);
        mapping.name = 'PCA';
    otherwise
        error('Unknown dimensionality reduction technique.');
end

end
