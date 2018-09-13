function Out = GillesF(I,radius,th)
% Input:     - I: A 2D image
%            - radius: radius of the region mask. Default: radius = 10
%            - th: only keep points above a threshold defined by th.
%              Default: th = 0.95
%
% Output:    - Out: A (1x6) vector containing 6 metrics calculated from
%                   detected Gilles points
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic 
% and Interventional Radiology, University Hospital of Tübingen, Germany 
% and the Institute of Signal Processing and System Theory University of 
% Stuttgart, Germany. Last modified: November 2016
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
% ************************************************************************

if ~exist('radius','var')
    radius = 10;
end    
if ~exist('th','var')
    th = 0.95;
end    

if any(~real(I))
   I = real(I); 
end   

%% extract Gilles points

% convert image to binary
BW = im2bw(I, graythresh(I));    
im = BW(:,:,1);

% define a region mask
mask = fspecial('disk',radius) > 0;

% compute local entropy
loc_entropy = entropyfilt(im,mask);

% find the local maxima
[~,~,local_max] = findLocalMaximum(loc_entropy,radius);

% keep only points above a threshold
[row,col] = find(local_max > th*max(local_max(:)));

% normalize coordinates for better comparison of different sized images
points_gilles = [row col]/numel(I);

%% feature extraction
points_gilles_num = size(points_gilles,1);
points_gilles_mean = mean(points_gilles);
points_gilles_std = std(points_gilles);
points_gilles_std2 = std2(points_gilles);



%% return features
Out = [points_gilles_num points_gilles_mean points_gilles_std points_gilles_std2];

end