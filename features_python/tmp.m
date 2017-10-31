I = csvread('testimg.csv');
typeflag.corr = false;
typeflag.texture = false;
typeflag.global = false;
typeflag.entropy = false;
addpath('ImFEATbox/GlobalFeatures/Geometrical/');
Out = GLCMF(I, typeflag);
csvwrite('matlab-out.csv', Out);