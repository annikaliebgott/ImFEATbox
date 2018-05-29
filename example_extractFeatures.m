%% Example script to carry out the feature extraction
%
% This script calls all feature extraction allgorithms provided by
% ImFEATbox.
%
% Output:   feat_vec: an [NxM] array of M features extracted from N 2D images.
%
% The basic steps to successfully carry out the feature extraction:
% 1.)   Import the images of which you want to extract the features.
%       Note: this script expects images to be cell arrays
% 2.)   Set typeflag to define which features you wish to extract.
%       For more information about typeflag: see README.txt
% 3.)   Set parameters needed for some feature extraction algorithms.
%       Note: default values are tuned for use with automated MR image
%       quality assessment, you might need to change them according to your
%       application
% 4.)   Choose wheather or not you wish to use additional tools
%       (preprocessing, visualization). 
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic
% and Interventional Radiology, University Hospital of Tuebingen, Germany
% and the Institute of Signal Processing and System Theory University of
% Stuttgart, Germany. Last modified: February 2017
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
%          thomas.kuestner@iss.uni-stuttgart.de
% ************************************************************************

%% add necessary paths
addpath(genpath([pwd,filesep,'features_matlab']));

%% import images to extract the features from
% Note: this script expects images to be saved as cell arrays
data = load(['dataset',filesep,'images.mat']);
images = data.images;


%% get available/implemented features
cFeatures = fExtractFeatures([],[], '-getFeatures');


%% extract all features 
feat_vector = [];
for iI=1:length(images)
    feat_vector = cat(1,fExtractFeatures(images{iI},[],'all'));
end
    
    
%% extract transformation features with feature parameters (specified in file)
sFeatureParameter = [pwd, filesep, 'parameters_ImFEATBox_def.m'];
feat_vector = [];
for iI=1:length(images)
    feat_vector = cat(1,fExtractFeatures(images{iI},[],'transform', sFeatureParameter));
end


%% extract specific features with feature parameters (specified in file)
sFeatureParameter = [pwd, filesep, 'parameters_ImFEATBox_def.m'];
feat_vector = [];
for iI=1:length(images)
    feat_vector = cat(1,fExtractFeatures(images{iI},[],{'intensity', 'mser'}, sFeatureParameter));
end

