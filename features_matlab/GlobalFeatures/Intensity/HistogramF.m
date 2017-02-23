
function Out = HistogramF(I,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.global: all features 
%                   + typeflag.texture: all features
%                   + typeflag.entropy: only features based on entropy
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x6) vector containing 6 metrics calculated from the
%                   image histogram
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

if ~exist('typeflag','var')
   typeflag.local = true; 
   typeflag.texture = true;
   typeflag.entropy = true;
end    

if typeflag.global || typeflag.texture
    typeflag.entropy = true;
end    

% Check for color image and convert to grayscale
if(numel(size(I))==3)
    if(size(I,3)==3)
        I = rgb2gray(I);
    end
end

% 256 gray scale Image
graylevels = (0:255); 

% Probability of occurence of gray values
Prob = histc(I(:),graylevels) ./ numel(I);

%% extract features

if (typeflag.texture || typeflag.global)
    % Histogram Mean _ Identifier code
    Mean = graylevels*Prob;
    
    % Histogram Standard Deviation
    Std = sqrt(((graylevels-Mean).^2)*Prob);
    
    % Histogram Skewness
    Skewness = 1/(Std^3)*((graylevels-Mean).^3)*Prob;
    
    % Histogram Kurtosis
    Kurtosis = 1/(Std^4)*(((graylevels-Mean).^4)*Prob-3);
    
    % Histogram Energy
    Energy = sum(Prob.^2);
end

% Histogram Entropy
% marginal entropies
Entropy = -sum(Prob(Prob ~= 0) .* log2(Prob(Prob ~= 0)));


%% Return feature vector

if (typeflag.texture || typeflag.global)
    Out = [Mean Std Skewness Kurtosis Energy Entropy];
else
    Out =  Entropy;
end

end
