I = csvread('testimg.csv');
typeflag.corr = false;
typeflag.texture = false;
typeflag.global = false;
typeflag.entropy = false;
addpath('ImFEATbox/GlobalFeatures/Intensity/');
Out = ImFEATbox/GlobalFeatures/Intensity/IntensityF(I, typeflag);
csvwrite('matlab-out.csv', Out);