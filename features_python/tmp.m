I = csvread('testimg.csv');
addpath('ImFEATbox/GlobalFeatures/Intensity/');
Out = SVDF(I);
csvwrite('matlab-out.csv', Out);