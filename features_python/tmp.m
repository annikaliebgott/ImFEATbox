I = csvread('testimg.csv');
typeflag.transform = false;
typeflag.corr = false;
typeflag.global = false;
addpath('ImFEATbox/GlobalFeatures/Transformation/');
Out = DCTF(I, typeflag);
csvwrite('matlab-out.csv', Out);