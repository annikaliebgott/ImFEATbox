function Out = SectorF(I)
% Input:     - I: A 2D image
%
%
% Output:    - Out: A (1x5) vector containing 5 metrics calculated from the
%              number of non-zero elements of different image sectors 
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

if any(~real(I))
   I = real(I); 
end   

%% break image into sectors

% convert image to binary image
bw = im2bw(I); 
[height, width] = size(bw); 

% get the center of the image
centroid = ceil([height, width]./2); 

% get indices of all pixels
[h,w] = meshgrid(1:height,1:width); 

% convert to polar coordinates relative to the image centre
[theta, ~] = cart2pol(h-centroid(1), w-centroid(1)); 

% calculate to which sector each pixel belongs
E = 0:4:180;
[~, SectorIdx] = histc(theta * (180/pi), E);


%% extract features

% count the number of non-zero elements for each sector
N = arrayfun(@(k) nnz(bw(SectorIdx==k)==0), 1:max(SectorIdx)); 

% calculate sum, mean, std, max and min of the number of non-zero pixels
% per sector
sum_N = sum(N);
mean_N = mean(N);
std_N = std(N);
max_N = max(N);
min_N = min(N);

%% return features
Out = [sum_N mean_N std_N max_N min_N];

end
