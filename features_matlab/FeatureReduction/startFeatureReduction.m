%% This script starts the feature reduction and saves the reduced feature matrix. 
% In case of unsupervised reduction algorithms, a feature matrix of
% size [N_samples x N_features] must be provided. To use the supervised
% reduction algorithms, a label vector of size [N_samples x 1]
% corresponding to the feature matrix must also be provided.
%
% The necessary parameters can be set below.
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
% Part of the implemented algorithms are based on or extracted from the 
% Matlab Toolbox for Dimensionality Reduction by Laurens van der Maaten
% (http://homepage.tudelft.nl/19j49), 
% (C) Laurens van der Maaten, Delft University of Technology
%
% Part of the implemented algorithms are based on or extracted from the
% MatlabFunc toolbox by Deng Cai
% (c) Deng Cai, Zhejiang University


%% set parameters
% choose which algorithm to use. Use EITHER 'all' OR 'unsupervised' OR
% 'supervised' or specify desired reduction techniques. 
% Implemented algorithms:
% 'PCA',  'DiffusionMaps', 'SPE', 'GPLVM', 'MDS', 'Sammon', 
% 'Isomap', 'LDA', 'IsoP', 'OLPP'

sAlgo = {'all'};
% 'all': use all reduction techniques
% 'unsupervised': use all unsupervised techniques  (PCA, DiffusionMaps, 
% SPE, GPLVM, MDS, Sammon, Isomap)
% 'supervised': use all supervised techniques (LDA, IsoP, 
% OLPP)

if strcmp(sAlgo, 'all')
   sAlgo = {'PCA',  'DiffusionMaps', 'SPE', 'GPLVM', 'MDS',...
       'Sammon', 'Isomap', 'LDA', 'IsoP', 'OLPP'};
end    
if strcmp(sAlgo, 'unsupervised')
   sAlgo = {'PCA',  'DiffusionMaps', 'SPE', 'GPLVM', 'MDS',...
       'Sammon', 'Isomap'}; 
end
if strcmp(sAlgo, 'supervised')
   sAlgo = {'LDA','IsoP', 'OLPP'} ;
end    

% choose how to handle NaN and Inf values. Options:
% 0: do nothing about NaN and Inf (causes problems with some algorithms).
% 1: to replace NaN and Inf with value specified in normvl.
% 2: replace only NaN with normvl.
% 3: replace only Inf with normvl.
ernan = 1;
normvl = 0;

% choose whether or not normalization should be performed
ifnorm = true;

% set maximum number of feature components that should be calculated by
% transformation algorithms
N_feat_max = 70;

% Variable fo Nearest Neigbour. Only needed for Isomap
NN = 12;


%% path settings
currpath = pwd;

% load reduction algorithms
addpath(genpath([currpath,filesep,'ReductionTechniques']));

% specify path and file name of the feature matrix you wish to reduce
sFilenameSamples = 'Samples.mat';
sFilenameLabels = 'Labels.mat';

% define path to save results
sSavePath = [cd,filesep,'SavedReductionResults'];

% define file name to save results
sSaveName = 'Results';

%% load data
% load feature matrix
data = load(sFilenameSamples);
varName = fieldnames(data);
Samples = getfield(data,varName{1});

if any(strcmp(sAlgo,'IsoP')) || any(strcmp(sAlgo,'OLPP')) ||...
        any(strcmp(sAlgo,'LDA'))
    try
        data = load(sFilenameLabels);
    catch
        disp('Error: file containing class labels could not be loaded.') 
        disp('Class labels need to be provided to perform supervised reduction techniques!')
        return
    end    
    varName = fieldnames(data);
    Labels = getfield(data,varName{1});
end

%% feature normalization
disp([num2str(sum(sum(isnan(Samples)))),' NaN Werte'])
disp([num2str(sum(sum(isinf(Samples)))),' Inf Werte'])

% handle NaN and Inf values
switch ernan 
    case 1
        Samples(isnan(Samples))=normvl;
        Samples(isinf(Samples))=normvl;
        disp('Nan and Inf values replaced by normvl')
    case 0
        disp('No Nan and Inf values replaced')
    case 2
        disp('No Inf values replaced, only NaN')
        Samples(isnan(Samples))=normvl;
    case 3
        disp('No Nan values replaced, only Inf')
        Samples(isinf(Samples))=normvl;
end

% normalize feature matrix using zscore
if ifnorm
    [Samples, ~, ~] = zscore(Samples);
end

%% feature reduction

for n = 1:length(sAlgo)
    
    Input.sAlgo = sAlgo{n};
    Input.N_feat_max = N_feat_max;
    Input.NN = NN;
    
    if any(strcmp(sAlgo,'IsoP')) || any(strcmp(sAlgo,'OLPP')) ||...
        any(strcmp(sAlgo,'LDA')) 
        Input.Labels = Labels;
    end
    
    % feature reduction
    Output = fFeatureReduction(Samples, Input);
    
    % save results
    save([sSavePath,filesep,sSaveName,sAlgo{n},'.mat'],'Output')
    
end
