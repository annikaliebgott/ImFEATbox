function Out = LoGF(I, sigma, N_blobs, typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.local: all features 
%                   + typeflag.texture: all features
%                   + typeflag.moments: only features based on moments
%              default: all features are being extracted
%              For more information see README.txt
%            - N_blobs: Desired number of points of interests. 
%              default: N_blobs = 120;
%            - sigma: a struct containing the fields begin, end and step.
%              default: sigma.begin = 2, sigma.end = 15, sigma.step = 1
%
%
% Output:    - Out: A (1x261) vector containing 261 metrics calculated from
%              points of interest based on Laplacian of Gaussian
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


% set default parameters
if ~exist('typeflag','var')
   typeflag.local = true; 
   typeflag.texture = true;
   typeflag.moments = true;
end 

if ~exist('N_blobs','var')
   N_blobs = 120; 
end    

if ~exist('sigma','var') || ~isfield('sigma','begin')
   sigma.begin = 2;
end    
if ~exist('sigma','var') || ~isfield('sigma','end')
   sigma.end = 15;
end    
if ~exist('sigma','var') || ~isfield('sigma','step')
   sigma.step = 1;
end    

% convert input image
I = double(I);

if any(~isreal(I))
   I = real(I); 
end   

% Laplacian of Gaussian parameters
sigma_array = sigma.begin:sigma.step:sigma.end;
sigma_nb = numel(sigma_array);

% variable
img_height = size(I,1);
img_width = size(I,2);

% calculate scale-normalized laplacian operator
snlo = zeros(img_height,img_width,sigma_nb);
for i=1:sigma_nb
    sigma = sigma_array(i);
    snlo(:,:,i) = sigma*sigma*imfilter(I,fspecial('log', floor(6*sigma+1), sigma),'replicate');
end

% search for local maxima
snlo_dil = imdilate(snlo,ones(3,3,3));
blob_candidate_index = find(snlo==snlo_dil);
blob_candidate_value = snlo(blob_candidate_index);
[~,index] = sort(blob_candidate_value,'descend');
blob_index = blob_candidate_index( index(1:min(N_blobs,numel(index))) );
[lig,col,sca] = ind2sub([img_height,img_width,sigma_nb],blob_index);
points = [lig,col,3*reshape(sigma_array(sca),[size(lig,1),1])];


%plot
%     imshow(img), hold on;
%     plot(points);


%% feature extraction

% calculate 2nd and 4th moments in x and y direction
m2_points_x = moment(points(:,1),2);
m2_points_y = moment(points(:,2),2);
m4_points_x = moment(points(:,1),4);
m4_points_y = moment(points(:,2),4);

sca_moment2 = moment(sca,2);
sca_moment4 = moment(sca,4);

if (typeflag.local || typeflag.texture)
    
    weight_points = sum(points(3));
    
    % mean distance between the points
    dist_points = zeros(1,length(points));
    for i = 1 : length(points)-1
        dist_points(i) = sqrt( (points(i+1,1) - points(i,1))^2 + (points(i+1,2) - points(i,2))^2);
    end
    dist_points_total = sum(dist_points);
    dist_points_mean = mean(dist_points);
    dist_points_min = min(dist_points);
    dist_points_max = max(dist_points);
    
    % relative location of points
    points_x = points(:,1)/size(I,1);
    points_y = points(:,2)/size(I,2);
    
    % center of gravity in x and y direction and for the points weighting
    av_x = mean2(points(:,1));
    av_y = mean2(points(:,2));
    av_w = mean2(points(:,3));
    
    % standard deviation in x and y direction and for the points weighting
    sv_x = std2(points(:,1));
    sv_y = std2(points(:,2));
    sv_w = std2(points(:,3));
    
    % feature extraction of blob candidate
    cand_mean = mean2(blob_candidate_value);
    cand_std = std2(blob_candidate_value);
    
    sca_mean = mean2(sca);
    sca_std = std2(sca);
end


%% return feature vector
if ~(typeflag.texture || typeflag.local)
    Out = [m2_points_x m2_points_y m4_points_x m4_points_y sca_moment4 sca_moment2];
else
    Out = [points_y.' points_x.' m2_points_x m2_points_y m4_points_x m4_points_y...
        weight_points dist_points_total dist_points_mean dist_points_min dist_points_max...
        av_x av_y av_w sv_x sv_y sv_w cand_mean cand_std...
        sca_mean sca_std sca_moment4 sca_moment2];
end

end
