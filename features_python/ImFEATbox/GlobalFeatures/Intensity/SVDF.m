function Out = SVDF(I)
% Input:     - I: A 2D image
%
%
% Output:    - Out: A (1x780) vector containing 780 metrics calculated
%              from singular value decomposition
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic
% and Interventional Radiology, University Hospital of Tuebingen, Germany
% and the Institute of Signal Processing and System Theory University of
% Stuttgart, Germany. Last modified: November 2016
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
% ************************************************************************


%% Calculate Singular Value Decomposition of the image

% convert image I to double
I = double(I);

% initialize feature variables
dia_elements = zeros(size(I,1),3);
eig_U = zeros(size(I,1),3);
eig_V = zeros(size(I,2),3);
det_U = zeros(1,3);
det_V = zeros(1,3);
trace_U = zeros(1,3);
trace_V = zeros(1,3);
rank_U = zeros(1,3);
rank_V = zeros(1,3);
median_eig_U = zeros(1,3);
median_eig_V = zeros(1,3);
max_eig_U = zeros(1,3);
max_eig_V = zeros(1,3);
mean_U = zeros(1,3);
mean_V = zeros(1,3);
mean_S = zeros(1,3);
std_U = zeros(1,3);
std_V = zeros(1,3);
std_S = zeros(1,3);
skewness_U = zeros(1,3);
skewness_V = zeros(1,3);
kurtosis_U = zeros(1,3);
kurtosis_V = zeros(1,3);

% Calculate the measures for 3 different orientations
for z=1 : 3

    if(z == 2)
        % rotate image by 90 degree
        I = imrotate(I,90,'bilinear','crop');
    elseif (z == 3)
        % rotate image by -90 degree
        I = imrotate(I,-180,'bilinear','crop');
    end

    % calculate singular value decomposition with diagonal matrix S and
    % unitary matrices U and V
    [U,S,V] = svd(I);

    %% feature extraction

    % calculate diagonal elements of matrix S
    for i = 1 : nnz(S)
        dia_elements(i,z) = S(i,i);
    end

    % eigen values of U and V
    eig_U(:,z) = eig(U);
    eig_V(:,z) = eig(V);
    
    % determinant of U and V
    det_U(z) = det(U);
    det_V(z) = det(V);

    % trace of U and V
    trace_U(z) = trace(U);
    trace_V(z) = trace(V);

    % rank of U and V
    rank_U(z) = rank(U);
    rank_V(z) = rank(V);

    % skewness of U and V
    skewness_U(z) = skewness(U(:));
    skewness_V(z) = skewness(V(:));

    % kurtosis of U and V
    kurtosis_U(z) = kurtosis(U(:));
    kurtosis_V(z) = kurtosis(V(:));

    % mean of U, V and S
    mean_U(z) = mean2(U);
    mean_V(z) = mean2(V);
    mean_S(z) = mean2(S);

    % standard deviation of U, V and S
    std_U(z) = std2(U);
    std_V(z) = std2(V);
    std_S(z) = std2(S);

    % median of eigen values of U and V
    median_eig_U(z) = median(eig_U(z,:));
    median_eig_V(z) = median(eig_V(z,:));

    % maximum of eigen values of U and V
    max_eig_U(z) = max(eig_U(z,:));
    max_eig_V(z) = max(eig_V(z,:));

end


%% return feature vector

Out = [reshape(dia_elements(1:40,:),1,numel(dia_elements(1:40,:)))...
    reshape(eig_U(1:100,:),1,numel(eig_U(1:100,:)))...
    reshape(eig_V(1:100,:),1,numel(eig_V(1:100,:)))...
    det_U det_V trace_U trace_V rank_U rank_V skewness_U skewness_V...
    kurtosis_U kurtosis_V mean_U mean_V mean_S std_U std_V std_S...
    median_eig_U median_eig_V max_eig_U max_eig_V];

end
