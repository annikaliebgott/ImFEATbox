I = csvread('testimg.csv');
InputParameters = {};
typeflag.corr = false;
typeflag.texture = false;
typeflag.global = true;
typeflag.entropy = false;
addpath('ImFEATbox/GlobalFeatures/Geometrical/');
Out = GLCMF(I, InputParameters, typeflag);
csvwrite('matlab-out.csv', Out);