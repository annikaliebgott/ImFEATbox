function Out = ConnectivityF(I)
% Input:     - I: A 2D image
%
%
% Output:    - Out: A (1x47) vector containing 47 metrics based on the 
%              connectity of detected objects in the image
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


if any(~isreal(I))
   I = real(I); 
end   

%% find connected components

% convert image to binary
BW = im2bw(I, graythresh(I));      

% calculate properties of regions (connected components) in the BW image
s = regionprops(BW,'area','centroid','Perimeter','Extrema','MajorAxisLength','MinorAxisLength');

% find connected components in binary image
CC = bwconncomp(BW);

% find weighted centroids of regions
s2 = regionprops(CC,I,'WeightedCentroid');


%% feature extraction

% Number of connected components (objects) in BW
N_co = CC.NumObjects;   

% Connectivity of the connected components
connectivity = CC.Connectivity;         
I_size = CC.ImageSize;

% area covered by objects in the binary image
total_area = bwarea(BW);     
covered_area = total_area / ( I_size(1)*I_size(2) );

% Euler number of binary image
eul = bweuler(BW);

% compactness, standardizad measure for circles
p=[s.Perimeter];
a=[s.Area];
compactness = p.^(2)./(4*pi*a);
compactness_mean = mean(compactness);
compactness_std = std(compactness);

% calculate mean and std of area, perimeter
a_mean = mean(a); 
a_std = std(a);
p_mean = mean(p);
p_std = std(p);

% Get centers and radii of the objects
centers = [s.Centroid];

% order center points column wise into x and y coordinates
centers = [centers(1:2:length(centers)-1)' centers(2:2:length(centers))']; 
centers_mean = mean(centers,1);
centers_std = std(centers,0,1);
major_axis_mean = mean([s.MajorAxisLength]);
minor_axis_mean = mean([s.MinorAxisLength]);
major_axis_std = std([s.MajorAxisLength]);
minor_axis_std = std([s.MinorAxisLength]);

% Get features from the weighted centroids
centers2 = [s2.WeightedCentroid];
centers2 = [centers2(1:2:length(centers2)-1)' centers2(2:2:length(centers2))'];
centers2_mean = mean(centers2,1);
centers2_std = std(centers2,0,1);


% Use features from extrema points in the region
% calculate mean value and the mean Euclidean distance to the object center
% for each of the 8 extrema points
extrema_sum = zeros(8,2);
extrema_dist = zeros(8,1);

for i = 1:length(a)
   extrema_sum(:,1) = extrema_sum(:,1) + s(i).Extrema(:,1);
   extrema_sum(:,2) = extrema_sum(:,2) + s(i).Extrema(:,2);
   extrema_dist(:) = extrema_dist + sum((s(i).Extrema - repmat(centers(i,:),8,1)).^2,2).^0.5;
end    

% mean value for each of the 8 extrema points
extrema_mean = extrema_sum./length(a);
extrema_mean = reshape(extrema_mean',16,1)';

% mean distance to object center for each of the 8 extrema points
extrema_dist_mean = (extrema_dist./length(a))';

Out = [N_co connectivity total_area covered_area eul... 
    compactness_mean compactness_std a_mean a_std p_mean p_std... 
    centers_mean centers_std centers2_mean centers2_std...
    major_axis_mean major_axis_std minor_axis_mean minor_axis_std...
    extrema_mean extrema_dist_mean];

end
