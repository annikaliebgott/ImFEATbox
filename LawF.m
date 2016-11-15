function Out = LawF(I)
% Input:     - I: A 2D image
%
%
% Output:    - Out: A (1x58) vector containing 58 metrics calculated from 
%                   the 2nd and 4th moments of the results of the image
%                   filtered by 25 different filters
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
% Law features are constructed from a set of five one-dimensional filters:
% edge, spot, ripple, wave and lowpass. By applying a one-dimensional 
% filter in the horizontal direction followed by a one-dimensional filter 
% in the vertical direction, this results in 25 different two-dimensional 
% filters.





%% filter image horizontal and vertical

%converte to double
I = double(I);

% define 5 filters
edge = [-1 -2 0 2 1];
low_pass = [1 4 6 4 1];
spots = [-1 0 2 0 1];
ripples = [1 -4 6 -4 1];
waves = [-1 2 0 -2 1];

% edge filter
ee = conv2(edge, edge, I, 'same');
le = conv2(edge, low_pass, I, 'same');
se = conv2(edge, spots, I, 'same');
re = conv2(edge, ripples, I, 'same');
we = conv2(edge, waves, I, 'same');

%spots
es = conv2(spots, edge, I, 'same');
ls = conv2(spots, low_pass, I, 'same');
ss = conv2(spots, spots, I, 'same');
rs = conv2(spots, ripples, I, 'same');
ws = conv2(spots, waves, I, 'same');

%ripples
er = conv2(ripples, edge, I, 'same');
lr = conv2(ripples, low_pass, I, 'same');
sr = conv2(ripples, spots, I, 'same');
rr = conv2(ripples, ripples, I, 'same');
wr = conv2(ripples, waves, I, 'same');

%waves
ew = conv2(waves, edge, I, 'same');
lw = conv2(waves, low_pass, I, 'same');
sw = conv2(waves, spots, I, 'same');
rw = conv2(waves, ripples, I, 'same');
ww = conv2(waves, waves, I, 'same');

%low pass
el = conv2(low_pass, edge, I, 'same');
ll = conv2(low_pass, low_pass, I, 'same');
sl = conv2(low_pass, spots, I, 'same');
rl = conv2(low_pass, ripples, I, 'same');
wl = conv2(low_pass, waves, I, 'same');



%% feature extraction

% determine the 2nd and 4th moment of the filtered images
m2 = [moment(ee(:),2) moment(se(:),2) moment(re(:),2) moment(we(:),2) moment(le(:),2)...
    moment(es(:),2) moment(ss(:),2) moment(rs(:),2) moment(ws(:),2) moment(ls(:),2)...
    moment(er(:),2) moment(sr(:),2) moment(rr(:),2) moment(wr(:),2) moment(lr(:),2)...
    moment(ew(:),2) moment(sw(:),2) moment(rw(:),2) moment(ww(:),2) moment(lw(:),2)...
    moment(el(:),2) moment(sl(:),2) moment(rl(:),2) moment(wl(:),2) moment(ll(:),2)];

m4 = [moment(ee(:),4) moment(se(:),4) moment(re(:),4) moment(we(:),4) moment(le(:),4)...
    moment(es(:),4) moment(ss(:),4) moment(rs(:),4) moment(ws(:),4) moment(ls(:),4)...
    moment(er(:),4) moment(sr(:),4) moment(rr(:),4) moment(wr(:),4) moment(lr(:),4)...
    moment(ew(:),4) moment(sw(:),4) moment(rw(:),4) moment(ww(:),4) moment(lw(:),4)...
    moment(el(:),4) moment(sl(:),4) moment(rl(:),4) moment(wl(:),4) moment(ll(:),4)];

% mean of the moments
mean_m2 = mean(m2);
mean_m4 = mean(m4);

% standard deviation of the moments
std_m2 = std2(m2);
std_m4 = std2(m4);

% maximum and minimum values of the moments
max_m2 = max(m2);
max_m4 = max(m4);
min_m2 = min(m2);
min_m4 = min(m4);


%% return feature vector

Out = [m2 m4 mean_m2 mean_m4 std_m2 std_m4 max_m2 max_m4 min_m2 min_m4];



end