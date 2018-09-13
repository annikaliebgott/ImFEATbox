function Out = GradientF(I,typeflag,gradtype)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features
%              based on desired characteristics:
%                   + typeflag.global: all features
%                   + typeflag.texture: all features
%                   + typeflag.gradient: all features
%                   + typeflag.entropy: only features based on entropy
%              default: all features are being extracted
%              For more information see README.txt
%            - gradtype: Struct of logicals to choose which type of
%              gradient should be applied:
%                   + gradtype.first: first order gradient
%                   + gradtype.second: second order gradient (Laplacian)
%              default: both types of gradients are used
%
%
% Output:    - Out: A (1x81) vector containing 81 metrics calculated from
%              image gradients
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
% ************************************************************************
%
% Implementation based on:  McGee et al.: "Image metric-based correction
%                           (Autocorrection) of motion effects: Analysis of
%                           image metrics", Journal of Magnetic Resonance
%                           Imaging, vol. 11, issue 2, p. 174-181, Feb 2000

%%
if ~exist('typeflag','var')
    typeflag.global = true;
    typeflag.texture = true;
    typeflag.entropy = true;
end

if ~exist('gradtype','var')
    gradtype.first = true;
    gradtype.second = true;
end

% Check for color image and convert to grayscale
if(numel(size(I))==3)
    if(size(I,3)==3)
        I = rgb2gray(I);
    end
end

Height = size(I,1);
Width = size(I,2);
n = Height * Width;

%% define gradients

if gradtype.first
    % gradient 1
    % gradient direction: x
    g1_x = [1; -1];
    % gradient direction: y
    g1_y = [1 -1];
    
    % gradient 2
    g2_x = [1; 0; -1];
    g2_y = [1 0 -1];
end

if gradtype.second
    % laplacian 1
    l1_x = [-1; 2; -1];
    l1_y = [-1 2 -1];
    
    % laplacian 2
    l2 = [-1 -2 -1; -2 12 -2; -1 -2 -1];
    
    % laplacian 3
    l3 = [0 -1 0; -1 4 -1; 0 -1 0];
    
    % laplacian 4
    l4 = [-1 -1 -1; -1 8 -1; -1 -1 -1];
end

%% convolve image
if gradtype.first
    G1_x = conv2(g1_x,I);
    G1_y = conv2(g1_y,I);
    G2_x = conv2(g2_x,I);
    G2_y = conv2(g2_y,I);
end

if gradtype.second
    L1_x = conv2(l1_x,I);
    L1_y = conv2(l1_y,I);
    L2 = conv2(l2,I);
    L3 = conv2(l3,I);
    L4 = conv2(l4,I);
end

%% extract features

% summed gradients
if gradtype.first
    G1_y_sum = sum(sum(abs(G1_y)));
    G1_x_sum = sum(sum(abs(G1_x)));
    G2_x_sum = sum(sum(abs(G2_x)));
    G2_y_sum = sum(sum(abs(G2_y)));
end
if gradtype.second
    L1_x_sum = sum(sum(abs(L1_x)));
    L1_y_sum = sum(sum(abs(L1_y)));
    L2_sum = sum(sum(abs(L2)));
    L3_sum = sum(sum(abs(L3)));
    L4_sum = sum(sum(abs(L4)));
end

% normalized gradients
if gradtype.first
    G1_x_norm = abs(G1_x)./G1_x_sum;
    G1_y_norm = abs(G1_y)./G1_y_sum;
    G2_x_norm = abs(G2_x)./G2_x_sum;
    G2_y_norm = abs(G2_y)./G2_y_sum;
end
if gradtype.second
    L1_x_norm = abs(L1_x)./L1_x_sum;
    L1_y_norm = abs(L1_y)./L1_y_sum;
    L2_norm = abs(L2)./L2_sum;
    L3_norm = abs(L3)./L3_sum;
    L4_norm = abs(L4)./L4_sum;
end

if typeflag.global || typeflag.texture || typeflag.gradient
    % sum of squared gradients
    if gradtype.first
        G1_x_2 = sum(sum((G1_x.^2)))/n;
        G1_y_2 = sum(sum((G1_y.^2)))/n;
        G2_x_2 = sum(sum((G2_x.^2)))/n;
        G2_y_2 = sum(sum((G2_y.^2)))/n;
    end
    if gradtype.second
        L1_x_2 = sum(sum((L1_x.^2)))/n;
        L1_y_2 = sum(sum((L1_y.^2)))/n;
        L2_2 = sum(sum((L2.^2)))/n;
        L3_2 = sum(sum((L3.^2)))/n;
        L4_2 = sum(sum((L4.^2)))/n;
    end
    
    % sum of 4th power of gradients
    if gradtype.first
        G1_x_4 = sum(sum((G1_x.^4)))/n;
        G1_y_4 = sum(sum((G1_y.^4)))/n;
        G2_x_4 = sum(sum((G2_x.^4)))/n;
        G2_y_4 = sum(sum((G2_y.^4)))/n;
    end
    if gradtype.second
        L1_x_4 = sum(sum((L1_x.^4)))/n;
        L1_y_4 = sum(sum((L1_y.^4)))/n;
        L2_4 = sum(sum((L2.^4)))/n;
        L3_4 = sum(sum((L3.^4)))/n;
        L4_4 = sum(sum((L4.^4)))/n;
    end
    
    % maximum of normalized gradients
    if gradtype.first
        G1_x_norm_max = max(G1_x_norm(:));
        G1_y_norm_max = max(G1_y_norm(:));
        G2_x_norm_max = max(G2_x_norm(:));
        G2_y_norm_max = max(G2_y_norm(:));
    end
    if gradtype.second
        L1_x_norm_max = max(L1_x_norm(:));
        L1_y_norm_max = max(L1_y_norm(:));
        L2_norm_max = max(L2_norm(:));
        L3_norm_max = max(L3_norm(:));
        L4_norm_max = max(L4_norm(:));
    end
    
    % standard deviation of normalized gradients
    if gradtype.first
        G1_x_norm_std = std2(G1_x_norm);
        G1_y_norm_std = std2(G1_y_norm);
        G2_x_norm_std = std2(G2_x_norm);
        G2_y_norm_std = std2(G2_y_norm);
    end
    if gradtype.second
        L1_x_norm_std = std2(L1_x_norm);
        L1_y_norm_std = std2(L1_y_norm);
        L2_norm_std = std2(L2_norm);
        L3_norm_std = std2(L3_norm);
        L4_norm_std = std2(L3_norm);
    end
    
    % mean of normalized gradients
    if gradtype.first
        G1_x_norm_mean = mean2(G1_x_norm);
        G1_y_norm_mean = mean2(G1_y_norm);
        G2_x_norm_mean = mean2(G2_x_norm);
        G2_y_norm_mean = mean2(G2_y_norm);
    end
    if gradtype.second
        L1_x_norm_mean = mean2(L1_x_norm);
        L1_y_norm_mean = mean2(L1_y_norm);
        L2_norm_mean = mean2(L2_norm);
        L3_norm_mean = mean2(L3_norm);
        L4_norm_mean = mean2(L3_norm);
    end
    
    % sum of normalized gradients squared
    if gradtype.first
        G1_x_norm_2 = sum(sum((G1_x_norm.^2)))/n;
        G1_y_norm_2 = sum(sum((G1_y_norm.^2)))/n;
        G2_x_norm_2 = sum(sum((G2_x_norm.^2)))/n;
        G2_y_norm_2 = sum(sum((G2_y_norm.^2)))/n;
    end
    if gradtype.second
        L1_x_norm_2 = sum(sum((L1_x_norm.^2)))/n;
        L1_y_norm_2 = sum(sum((L1_y_norm.^2)))/n;
        L2_norm_2 = sum(sum((L2_norm.^2)))/n;
        L3_norm_2 = sum(sum((L3_norm.^2)))/n;
        L4_norm_2 = sum(sum((L4_norm.^2)))/n;
    end
    
    % sum of normalized gradients to 4th power
    if gradtype.first
        G1_x_norm_4 = sum(sum((G1_x_norm.^4)))/n;
        G1_y_norm_4 = sum(sum((G1_y_norm.^4)))/n;
        G2_x_norm_4 = sum(sum((G2_x_norm.^4)))/n;
        G2_y_norm_4 = sum(sum((G2_y_norm.^4)))/n;
    end
    if gradtype.second
        L1_x_norm_4 = sum(sum((L1_x_norm.^4)))/n;
        L1_y_norm_4 = sum(sum((L1_y_norm.^4)))/n;
        L2_norm_4 = sum(sum((L2_norm.^4)))/n;
        L3_norm_4 = sum(sum((L3_norm.^4)))/n;
        L4_norm_4 = sum(sum((L4_norm.^4)))/n;
    end
end

% marginal entropies of normalized gradients
if gradtype.first
    G1_x_E = -sum(G1_x_norm(G1_x_norm ~= 0) .* log2(G1_x_norm(G1_x_norm ~= 0)));
    G1_y_E = -sum(G1_y_norm(G1_y_norm ~= 0) .* log2(G1_y_norm(G1_y_norm ~= 0)));
    G2_x_E = -sum(G2_x_norm(G2_x_norm ~= 0) .* log2(G2_x_norm(G2_x_norm ~= 0)));
    G2_y_E = -sum(G2_y_norm(G2_y_norm ~= 0) .* log2(G2_y_norm(G2_y_norm ~= 0)));
end
if gradtype.second
    L1_x_E = -sum(L1_x_norm(L1_x_norm ~= 0) .* log2(L1_x_norm(L1_x_norm ~= 0)));
    L1_y_E = -sum(L1_y_norm(L1_y_norm ~= 0) .* log2(L1_y_norm(L1_y_norm ~= 0)));
    L2_E = -sum(L2_norm(L2_norm ~= 0) .* log2(L2_norm(L2_norm ~= 0)));
    L3_E = -sum(L3_norm(L3_norm ~= 0) .* log2(L3_norm(L3_norm ~= 0)));
    L4_E = -sum(L4_norm(L4_norm ~= 0) .* log2(L4_norm(L4_norm ~= 0)));
end


%% Build output vector Out
if typeflag.global || typeflag.texture || typeflag.gradient
    if gradtype.first && gradtype.second
        Out = [G1_x_sum G1_y_sum G2_x_sum G2_y_sum L1_x_sum L1_y_sum L2_sum L3_sum L4_sum...
            G1_x_2 G1_y_2 G2_x_2 G2_y_2 L1_x_2 L1_y_2 L2_2 L3_2 L4_2...
            G1_x_4 G1_y_4 G2_x_4 G2_y_4 L1_x_4 L1_y_4 L2_4 L3_4 L4_4...
            G1_x_norm_max G1_y_norm_max G2_x_norm_max G2_y_norm_max...
            L1_x_norm_max L1_y_norm_max L2_norm_max L3_norm_max L4_norm_max...
            G1_x_norm_std G1_y_norm_std G2_x_norm_std G2_y_norm_std...
            L1_x_norm_std L1_y_norm_std L2_norm_std L3_norm_std L4_norm_std...
            G1_x_norm_mean G1_y_norm_mean G2_x_norm_mean G2_y_norm_mean...
            L1_x_norm_mean L1_y_norm_mean L2_norm_mean L3_norm_mean L4_norm_mean...
            G1_x_norm_2 G1_y_norm_2 G2_x_norm_2 G2_y_norm_2...
            L1_x_norm_2 L1_y_norm_2 L2_norm_2 L3_norm_2 L4_norm_2...
            G1_x_norm_4 G1_y_norm_4 G2_x_norm_4 G2_y_norm_4...
            L1_x_norm_4 L1_y_norm_4 L2_norm_4 L3_norm_4 L4_norm_4...
            G1_x_E G1_y_E G2_x_E G2_y_E L1_x_E L1_y_E L2_E L3_E L4_E];
    else if gradtype.first
            Out = [G1_x_sum G1_y_sum G2_x_sum G2_y_sum...
            G1_x_2 G1_y_2 G2_x_2 G2_y_2...
            G1_x_4 G1_y_4 G2_x_4 G2_y_4...
            G1_x_norm_max G1_y_norm_max G2_x_norm_max G2_y_norm_max...
            G1_x_norm_std G1_y_norm_std G2_x_norm_std G2_y_norm_std...
            G1_x_norm_mean G1_y_norm_mean G2_x_norm_mean G2_y_norm_mean...
            G1_x_norm_2 G1_y_norm_2 G2_x_norm_2 G2_y_norm_2...
            G1_x_norm_4 G1_y_norm_4 G2_x_norm_4 G2_y_norm_4...
            G1_x_E G1_y_E G2_x_E G2_y_E];
        else
            Out = [L1_x_sum L1_y_sum L2_sum L3_sum L4_sum...
            L1_x_2 L1_y_2 L2_2 L3_2 L4_2...
            L1_x_4 L1_y_4 L2_4 L3_4 L4_4...
            L1_x_norm_max L1_y_norm_max L2_norm_max L3_norm_max L4_norm_max...
            L1_x_norm_std L1_y_norm_std L2_norm_std L3_norm_std L4_norm_std...
            L1_x_norm_mean L1_y_norm_mean L2_norm_mean L3_norm_mean L4_norm_mean...
            L1_x_norm_2 L1_y_norm_2 L2_norm_2 L3_norm_2 L4_norm_2...
            L1_x_norm_4 L1_y_norm_4 L2_norm_4 L3_norm_4 L4_norm_4...
            L1_x_E L1_y_E L2_E L3_E L4_E];
        end
    end
else if gradtype.first && gradtype.second
    Out = [G1_x_E G1_y_E G2_x_E G2_y_E L1_x_E L1_y_E L2_E L3_E L4_E];
    else if gradtype.first
            Out = [G1_x_E G1_y_E G2_x_E G2_y_E];
        else
            Out = [L1_x_E L1_y_E L2_E L3_E L4_E];
        end
    end    
end
