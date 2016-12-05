function Out = LineProfileF(I,plotflag,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.local: all features 
%                   + typeflag.texture: all features
%                   + typeflag.corr: only features based on correlation
%                   + typeflag.moments: only features based on moments
%              default: all features are being extracted
%              For more information see README.txt
%            - plotflag: logical flag to enable/disable visualization.
%              default: plotflag = false;
%
%
% Output:    - Out: A (1x122) vector containing 122 metrics calculated from
%                   an intensity profile of the image
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
   typeflag.corr = true;
   typeflag.moments = true;
end    
if~exist('plotflag','var')
    plotflag = false;
end

%% Extract the intensity profile along different line segments

x_max = double(size(I,2));
y_max = double(size(I,1));
x_half = double((size(I,2) / 2));
y_half = double((size(I,1) / 2));

% diagonal lines with fitting size for the input image
x1 = [1 x_max];
y1 = [1 y_max];
x2 = [x_max 1];
y2 = [1 y_max];

% vertical line through the center of the image
x3 = [x_half x_half];
y3 = [1 y_max];

% horizontal line through the center of the image
x4 = [1 x_max];
y4 = [y_half y_half];

% view line profile
if plotflag
    x = [x1 x2];
    y = [y1 y2];
    improfile(I,x,y),grid on;
end

% extract intensity profiles along the defined lines
c1 = improfile(I,x1,y1);
c2 = improfile(I,x2,y2);
c3 = improfile(I,x3,y3);
c4 = improfile(I,x4,y4);

l1= size(c1,1);
l2= size(c2,1);
l3= size(c3,1);
l4= size(c4,1);

% Use zero-padding to get the same size for all arrays
m = max([l1 l2 l3 l4]);
c1 = padarray(c1, (m-l1), 'post');
c2 = padarray(c2, (m-l2), 'post');
c3 = padarray(c3, (m-l3), 'post');
c4 = padarray(c4, (m-l4), 'post');


%% Feature extraction

% linear correlation
if (typeflag.corr || typeflag.local || typeflag.texture)
    roh12 = corr(c1,c2);
    roh13 = corr(c1,c3);
    roh14 = corr(c1,c4);
    roh23 = corr(c2,c3);
    roh24 = corr(c4,c2);
    roh34 = corr(c3,c4);
end

% central moments
if (typeflag.moments || typeflag.local || typeflag.texture)
    m1 = moment(c1,2);
    m2 = moment(c2,2);
    m3 = moment(c3,2);
    m4 = moment(c4,2);
    
    m1_5 = moment(c1,5);
    m2_5 = moment(c2,5);
    m3_5 = moment(c3,5);
    m4_5 = moment(c4,5);
end

if (typeflag.local || typeflag.texture)
    % average value of array
    mean1 = mean(c1);
    mean2 = mean(c2);
    mean3 = mean(c3);
    mean4 = mean(c4);
    
    % standard deviation
    sd1 = std(c1);
    sd2 = std(c2);
    sd3 = std(c3);
    sd4 = std(c4);
    
    % Percentiles of a data set
    pr1 = prctile(c1,48);   
    pr2 = prctile(c2,48);
    pr3 = prctile(c3,48);
    pr4 = prctile(c4,48);
    
    % fast fourier transform
    fft1 = fft(c1);
    fft2 = fft(c2);
    fft3 = fft(c3);
    fft4 = fft(c4);
    
    % power of spectrum
    p1 = abs(sum(power(fft1,2)) / length(fft1));
    p2 = abs(sum(power(fft2,2)) / length(fft2));
    p3 = abs(sum(power(fft3,2)) / length(fft3));
    p4 = abs(sum(power(fft4,2)) / length(fft4));
    
    % derivates
    dc1 = diff(c1);
    dc2 = diff(c2);
    dc3 = diff(c3);
    dc4 = diff(c4);
    
    % average value of array of the derivates
    mean_dc1 = mean(dc1);
    mean_dc2 = mean(dc2);
    mean_dc3 = mean(dc3);
    mean_dc4 = mean(dc4);
    
    % standard deviation of the derivates
    sd_dc1 = std(dc1);
    sd_dc2 = std(dc2);
    sd_dc3 = std(dc3);
    sd_dc4 = std(dc4);
    
    % fast fourier transform of the derivates
    fft_dc1 = fft(dc1);
    fft_dc2 = fft(dc2);
    fft_dc3 = fft(dc3);
    fft_dc4 = fft(dc4);
    
    % power of spectrum of the derivates
    p_dc1 = abs(sum(power(fft_dc1,2)) / length(fft_dc1));
    p_dc2 = abs(sum(power(fft_dc2,2)) / length(fft_dc2));
    p_dc3 = abs(sum(power(fft_dc3,2)) / length(fft_dc3));
    p_dc4 = abs(sum(power(fft_dc4,2)) / length(fft_dc4));
    
end


%% Return feature vector

if ~(typeflag.local || typeflag.texture)
    if ~typeflag.moments
        Out = [roh12 roh14 roh24 roh13 roh23 roh34];
    elseif ~typeflag.corr
        Out = [m1 m2 m3 m4 m1_5 m2_5 m3_5 m4_5];
    else
        Out = [roh12 roh14 roh24 roh13 roh23 roh34...
            m1 m2 m3 m4 m1_5 m2_5 m3_5 m4_5];
    end
else
    Out = [roh12 roh14 roh24 roh13 roh23 roh34...
        mean1 mean2 mean3 mean4...
        sd1 sd2 sd3 sd4...
        pr1  pr2  pr3 pr4...
        m1 m2 m3 m4 m1_5 m2_5 m3_5 m4_5...
        p1 p2 p3 p4...
        fft1(1:10,1).' fft2(1:10,1).' fft3(1:10,1).' fft4(1:10,1).'...
        mean_dc1 mean_dc2 mean_dc3 mean_dc4...
        sd_dc1 sd_dc2 sd_dc3 sd_dc4...
        p_dc2 p_dc1 p_dc3 p_dc4...
        fft_dc1(1:10,1).' fft_dc2(1:10,1).' fft_dc3(1:10,1).' fft_dc4(1:10,1).'];
end

end