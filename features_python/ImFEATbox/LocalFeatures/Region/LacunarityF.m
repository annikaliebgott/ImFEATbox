function Out = LacunarityF(I,typeflag,l_min,l_max, stepsize)
% Input:     - I: A 2D real image
%            - typeflag: Struct of logicals to permit extracting features
%              based on desired characteristics:
%                   + typeflag.local: all features
%                   + typeflag.texture: all features
%                   + typeflag.corr: only features based on correlation
%              default: all features are being extracted
%              For more information see README.txt
%            - l_min: A scalar containing the minimum length of the box
%              sliding over the image. Default: l_min = 10
%            - l_max: A scalar indicating the maximum length of the box
%              sliding over the image. Default: l_max = 30
%            - stepsize: A scalar indicating the desired step size between
%              different box sizes. Default: stepsize = 1
%
%
% Output:    - Out: A 1x6 vector
%              containing metrics calculated from the image lacunarity
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

if ~exist('typeflag','var')
    typeflag.local = true;
    typeflag.texture = true;
    typeflag.corr = true;
end

if ~exist('l_min','var')
    l_min = 10;
end

if ~exist('l_max','var')
    l_max = 30;
end

if ~exist('stepsize','var')
    stepsize = 1;
end


%% catch errors

while any(size(I) < l_min)
    disp('Warning: l_min is too large for your image.')
    size_I_min = num2str(min(size(I)));
    l_min = input(['Enter a new value for l_min smaller than ',size_I_min,': ']);
end    
while l_min > l_max
    disp('Warning: l_max must be larger than l_min.')  
    l_min_str = num2str(l_min);
    l_max = input(['Enter new value for l_max larger than ',l_min_str,': ']);
end    
if any(size(I) < l_max)
    l_max = min(size(I));
    disp(['Warning: l_max is larger than one or both sides of your image. l_max has been set to min(size(I)) = ',num2str(l_max),'.']);    
end
while stepsize > l_max - l_min
    disp('Warning: stepsize is to big for your current values of l_min and l_max.');
    stepsize_max = num2str(l_max-l_min);
    stepsize = input(['Enter a new stepsize smaller or equal to ',stepsize_max,': ']);
end

if(~isreal(I))
%     warning('LacunarityF(): Only processing of real input data');
    I = real(I);
end

%% calculate lacunarity
% convert image I to binary image
    I =  im2bw(double(I));
    
    l_vec = (l_min:stepsize:l_max);
    
    % create lacunarity array
    lacunarity = zeros(1,length(l_vec));
    
    % calculate lacunarity for each desired box length
    for i = 1:length(l_vec)
        % calculate the number of zero elements (= black pixels) by 
        % employing a 2 dimensional convolution on the image
        temp = conv2(double(I==0), ones(l_vec(i)), 'valid');
       
        % calculate the average value and standard derivation of temp
        mean_temp = mean2(temp);
        std_temp = std2(temp);
        
        % calculate lacunarity
        lacunarity(i) = (std_temp/mean_temp)^2;
    end
    
    %% calculate features from lacunarity
    
    if typeflag.local || typeflag.texture
        % the more uniform the texture is, the smaller the mean
        lacunarity_max = max(lacunarity);
        lacunarity_min = min(lacunarity);
        lacunarity_mean = mean(lacunarity);
        lacunarity_std = std(lacunarity);
    end
    
    % calculate cross-covariance
    lacunarity_xcov = xcov(lacunarity);
    lacunarity_xcov_mean = mean(lacunarity_xcov);
    lacunarity_xcov_std = std(lacunarity_xcov);
    
    %% return feature vector
    
    if typeflag.local || typeflag.texture
        Out = [lacunarity_min lacunarity_max lacunarity_std lacunarity_mean lacunarity_xcov_std lacunarity_xcov_mean];
    else
        Out = [lacunarity_xcov_std lacunarity_xcov_mean];
    end

end
