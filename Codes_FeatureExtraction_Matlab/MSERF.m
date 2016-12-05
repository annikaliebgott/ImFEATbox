function Out = MSERF(I,plotflag)
% Input:     - I: A 2D image
%            - plotflag: logical flag to enable/disable visualization
%
% Output:    - Out: A (1x15) vector containing 15 metrics based on detected
%              maximally stable extremal regions in the image
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
% Reference:   	J. Matas, O. Chum, M. Urban and T. Pajdla. Robust Wide 
%               Baseline Stereo from Maximally Stable Extremal Regions. 
%               In David Marshall and Paul L. Rosin, editors, Proceedings 
%               of the British Machine Conference, pages 36.1-36.10. 
%               BMVA Press, September 2002. 


%% Detection of MSERs

if ~exist('plotflag','var')
   plotflag = false; 
end    

% convert image
I = uint8(I); 

% detect and store regions
regions = detectMSERFeatures(I);

%% Feature extraction

% total number of regions
num = regions.Count;

if regions.length > 0
    
    % determine center of gravity for all points on x- and y-axis, mean
    % orientation and mean area covered for all regions, mean coordinates 
    % of MSERs, std of coordinates of MSERs, sum of x- and y-coordinates
    a = regions.Location;
    o = regions.Orientation;
    axes = regions.Axes;
    total_area = size(I,1)*size(I,2);
    
    % center of gravity for both axes
    x_gravity_allpoints = sum(a(:,1)/regions.length);
    y_gravity_allpoints = sum(a(:,2)/regions.length);
    
    % mean orientation
    orientation_mean = sum(o)/regions.length;
    
    % mean area covered
    area_covered = sum(pi*axes(:,1).*axes(:,2)/total_area);
    
    % mean coordinates of MSERs
    axes_mean = mean2(axes);
    
    % std of coordinates of MSERs
    axes_std = std2(axes);
    
    % sum of coordinates of MSERs
    axes_sum_x = sum(axes(:,1));
    axes_sum_y = sum(axes(:,2));
    
    % features from pixel list: length, mean of x- and 
    % y-coordinates, std of x- and y-coordinates, difference of std of y-
    % and y-coordinates
    pixlist = cell2mat(regions.PixelList);  
    pixlist = double(pixlist);
    pixlist_length = length(pixlist);
    pixlist_mean_x = mean(pixlist(:,1));
    pixlist_mean_y = mean(pixlist(:,2));
    pixlist_std_x = std(pixlist(:,1));
    pixlist_std_y = std(pixlist(:,2));
    pixlist_std_diff = std(double(pixlist(:,1))) - std(double(pixlist(:,2)));
    
else
    
    x_gravity_allpoints = 0;
    y_gravity_allpoints = 0;
    orientation_mean = 0;
    area_covered = 0;
    axes_mean = 0;
    axes_std = 0;
    axes_sum_x = 0;
    axes_sum_y = 0;
    pixlist_length = 0;
    pixlist_mean_x = 0;
    pixlist_mean_y = 0;
    pixlist_std_x = 0;
    pixlist_std_y = 0;
    pixlist_std_diff =0;
 
end    

% display the centroids and axes of detected regions
if plotflag
    imshow(I); hold on;
    plot(regions); hold all;
    plot(x_gravity_allpoints, y_gravity_allpoints,'*r');
    %mark selectred regions colored
    plot(regions, 'showPixelList', true);     
end

%% return feature vector
Out = [num x_gravity_allpoints y_gravity_allpoints area_covered...
    orientation_mean axes_mean axes_std axes_sum_x axes_sum_y...
    pixlist_length pixlist_mean_x pixlist_mean_y...
    pixlist_std_x pixlist_std_y pixlist_std_diff];

end      


