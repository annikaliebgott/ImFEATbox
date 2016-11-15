function Out = SalientRegionF(I,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.texture: all features
%                   + typeflag.moments: only features based on moments
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x8) vector containing 8 metrics calculated from
%              detected salient regions
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
% Implementation based on:  R. Achanta, F. Estrada, P. Wils and S. 
%                           Süsstrunk, Salient Region Detection and 
%                           Segmentation, International Conference on 
%                           Computer Vision Systems (ICVS '08), Vol. 5008, 
%                           Springer Lecture Notes in Computer Science, 
%                            pp. 66-75, 2008.


if ~exist('typeflag','var')
   typeflag.texture = true;
   typeflag.moments = true;
end    


%% saliency calculation

x_length = size(I,2);
y_length = size(I,1);

if x_length > y_length
    shortest_side = y_length;
else
    shortest_side = x_length;
end

% Preallocate final saliency map M, contrast based saliency value C_help
M = zeros(y_length, x_length);
C = zeros(y_length, x_length);

% Determine haar feature vector for every pixel. Border areas are not
% considered.
C_horizontal = zeros(y_length, x_length);
C_vertical = zeros(y_length, x_length);

% set border
b = 14;
b2 = 12;
% [horz, vert] = Haarfeature_new(I,j,i,length), length = 1 (1 pixel). For
% other lengths modification of HaarF.m is necessary
for i = b : (y_length - b)
    for j = b : (x_length - b)
        [horz, vert] = HaarF(I,j,i,1,b2);
        C_horizontal(i,j) = horz;
        C_vertical(i,j) = vert;
    end
end

% various scales for the detection filter
% w = [2 4 8]: higher or smaller values would result in bad saliency values
for w = [2 4 8]
    
    % determine the width of the detection filter
    w_R2 = round(shortest_side / w);
    
    % move the detection filter in a raster scan fashion over the image
    for i = (w_R2) : (y_length-w_R2)
        for j = (w_R2) : (x_length-w_R2)
            
            % inner region R1 is usually chosen to be one pixel,
            % if noisy  choose NxN pixles
            R1 = 1;
            
            % Haar features
            % feature vector for R1
            v_q = sum([C_horizontal(i,j); C_vertical(i,j)]);
            % feature vector for R2
            w_h = floor(w_R2/2);
            v_p = sum([C_horizontal((i-w_h):(i+w_h),(j-w_h):(j+w_h));...
                C_vertical((i-w_h):(i+w_h),(j-w_h):(j+w_h))]);
            
            % correction term
            v_p = v_p - v_q;
            
            % determine distance between the average vectors of pixel
            % features of the inner region R1 and outer region R2
            % euclidean distance weight functions
            C(i,j) = sqrt(sum((v_p/w_R2 - v_q/R1).^2));
            
        end
    end
    
    %final saliency map is determined as sum of saliency values across the
    %scale
    M = double(M) + double(C);
    
end

if typeflag.texture
    % search for maxima in the final saliency map M
    max_position = [1,1];
    for i = 1 : y_length
        for j = 1 : x_length
            if (M(i,j) > M(max_position(1,1),max_position(1,2)))
                max_position = [i,j];
            end
        end
    end
    
    % relative max position
    max_position = max_position/(x_length*y_length);
    
    [M_max_column , ~] = max(M);
    [M_max_row , ~] = max(M,[],2);
    
    % average value
    mean_colum = mean(M_max_column);
    mean_row = mean(M_max_row);
    
    % standard deviation
    st_colum = std(M_max_column);
    st_row = std(M_max_row);
end

% moments
mom_M_2 = mean(moment(M,2));
mom_M_4 = mean(moment(M,4));

%     imshow(I), hold on;
%     plot(I,'+','red');

%% return features
if typeflag.texture
    Out = [mean_colum mean_row st_colum st_row max_position mom_M_2 mom_M_4];
else
    Out = [mom_M_2 mom_M_4];
end
end