function Out = SkeletonizationF(I,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.global: all features 
%                   + typeflag.transform: all features
%                   + typeflag.corr: only features based on correlation
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x17) vector containing 17 metrics based on
%              skeletonization transformation
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
   typeflag.global = true; 
   typeflag.transform = true;
   typeflag.corr = true;
end    

%% Transformation
% converte image

% im2bw can't process complex input values
I = double(real(I)); 

BW1 = im2bw(double(real(I)));

% skeletonization
BW2 = bwmorph(BW1,'skel',Inf);

% find the perimeter pixels in binary image
BW3 = bwperim(BW1);



%% feature extraction

if (typeflag.global || typeflag.transform)
% The density of a sparse matrix
dens_skelet = nnz(BW2)/numel(BW2);
length_line_skelet = nnz(BW2);
dens_peri = nnz(BW3)/numel(BW3);
length_line_peri = nnz(BW3);

% center of gravity
[x2, y2] = find(BW2);
c_g2 = [mean(x2) mean(y2)];
[x3, y3] = find(BW3);
c_g3 = [mean(x3) mean(y3)];

% mean value and std
diff = BW2 - BW3;
mean_2 = mean2(BW2);
mean_3 = mean2(BW3);
mean_d = abs(mean2(diff));
std_2 = std2(BW2);
std_3 = std2(BW3);
std_d = abs(std2(diff));
end

% 2D cross correlation
corr_coef = corr2(BW2,BW3);
corr = xcorr2(double(BW2),double(BW3));
corr_mean = mean2(corr);
corr_std = std2(corr);


%% return feature vector
if (typeflag.global || typeflag.transform)
Out = [dens_skelet length_line_skelet dens_peri length_line_peri c_g2 c_g3...
    mean_2 mean_3 mean_d std_2 std_3 std_d...
    corr_coef corr_mean corr_std];
else
   Out = [corr_coef corr_mean corr_std];
end    

end
