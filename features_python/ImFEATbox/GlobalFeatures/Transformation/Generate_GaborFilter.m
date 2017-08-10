function GaborResult = Generate_GaborFilter(scale,orientation,I, FlagShow)
% This function is called by GaborFilterF.m. It generates different gabor 
% filters with different scales and different orientations. Its output is a 
% (scale x orientation) cell whose elements are Gabor filtered images.
%
% Input:     - I: A 2D image
%            - scale: Number of scales to create the filters
%            - orientation: Number of orientations to create the filters
%            - Flagshow: A logical flag to enable/disable visualization.
%                        default: false
%
%
% Output:    - GaborResult: A (scale x orientation) cell so that each 
%                           element includes a matrix of 40x40 representing
%                           the corresponding filter
%
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
% Implementation based on:  J. K. Kamarainen: "Gabor features in image 
%                           analysis", 3rd International Conference on 
%                           Image Processing Theory, Tools and Applications
%                           (IPTA), Istanbul, 2012, pp. 13-14.


%%
% Create scale*orientation gabor filters
GaborFilter = cell(scale,orientation);

for i = 1:scale
    
    f_max = (0.25)/((sqrt(2))^(i-1));
    alpha = f_max/(sqrt(2));
    beta = f_max/sqrt(2);
    
    for j = 1:orientation
        theta = ((j-1)/orientation)*(2*pi);
        filter = zeros(39,39);
        y = 1:39;
        for x = 1:39
           xprime = (x-20)*cos(theta)+(y-20).*sin(theta);
           yprime = -(x-20)*sin(theta)+(y-20).*cos(theta);
           filter(x,:) = ((f_max)^2/(pi*(sqrt(2))*(sqrt(2))))*...
               exp(-((alpha^2).*(xprime.^2)+(beta^2).*(yprime.^2))).*exp(1i*2*pi*f_max.*xprime);
        end
        GaborFilter{i,j} = filter;
    end
end

if FlagShow
    %% Show Gabor filters
    
    % Show magnitudes of Gabor filters:
    figure('Name','Magnitudes of Gabor filters');
    for i = 1:scale
        for j = 1:orientation
            subplot(scale,orientation,(i-1)*orientation+j);
            imshow(abs(GaborFilter{i,j}),[]);
        end
    end
    
    % Show real parts of Gabor filters:
    figure('Name','Real parts of Gabor filters');
    for i = 1:scale
        for j = 1:orientation
            subplot(scale,orientation,(i-1)*orientation+j);
            imshow(real(GaborFilter{i,j}),[]);
        end
    end
end


% Convolve the Gabor filters with the given image
GaborResult = cell(scale,orientation);
for i = 1:scale
    for j = 1:orientation
        GaborResult{i,j} = conv2(I,GaborFilter{i,j},'same');
    end
end

