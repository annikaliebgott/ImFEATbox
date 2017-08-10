function Out = IntensityF(I,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.global: all features 
%                   + typeflag.texture: all features 
%                   + typeflag.corr: only features based on correlation
%                   + typeflag.entropy: only features based on entropy
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x7) vector containing 7 metrics based on intensity 
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
%
% Implementation based on the intensity-based features of the paper by 
% McGee et al.: "Image metric-based correction (Autocorrection) of motion 
% effects: Analysis of image metrics" (Journal of Magnetic Resonance
% Imaging, vol. 11, no. 2, pp. 174–181, 2000.)

if ~exist('typeflag','var')
    typeflag.global = true;
    typeflag.texture = true;
    typeflag.corr = true;
    typeflag.entropy = true;
end    

if typeflag.global || typeflag.texture
    typeflag.corr = true;
    typeflag.entropy = true;
end    

% Check for color image and convert to grayscale image
if(numel(size(I))==3)
    if(size(I,3)==3)
        I = rgb2gray(I);
    end
end

Height = size(I,1);
Width = size(I,2);
n = Height * Width;

%% Extract features

if (typeflag.global || typeflag.texture || typeflag.corr)
    % AutoCorrelation 1
    p2 = zeros(Height,Width-1);
    p1 = sum(I(:).^2);
    k = (1:Width-1);
    for i = 1:Height
        p2(i,k) = I(i,k).*I(i,k+1);
    end
    ACORR = p1 - sum(p2(:));
    
    % AutoCorrelation 2
    p3 = zeros(Height,Width-2);
    l = (1:Width-2);
    for i = 1:Height
        p3(i,l) = I(i,l).*I(i,l+2);
    end
    ACORR2 = sum(p2(:)) - sum(p3(:));
end

if (typeflag.global || typeflag.texture || typeflag.entropy)
    % Entropy
    H1 = double(I) ./ (sum((double(I(:)).^2)));
    % marginal entropies
    E = -sum(H1(H1 ~= 0) .* log(H1(H1 ~= 0)));
end

if (typeflag.global || typeflag.texture)
    
    %Standard Deviation
    STD = std(I(:));
    
    % Cube of Normalized Intensities
    NI3 = sum((I(:)./sum(I(:))).^3);
    
    % 4th power of Normalized Intensities
    NI4 = sum((I(:)./sum(I(:))).^4);
    
    % Squared Intensities
    I2 = sum((I(:)./n).^2);
end

%% return feature vector
if ~(typeflag.texture || typeflag.global)
    if ~typeflag.corr
        Out = E;
    elseif ~typeflag.entropy
        Out = [ACORR ACORR2];
    else
        Out = [ACORR ACORR2 E];
    end
else
    Out = [STD ACORR ACORR2 E NI3 NI4 I2];
end


end