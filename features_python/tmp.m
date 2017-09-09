I = csvread('testimg.csv');
gradtype.first = false;
gradtype.second = false;
typeflag.texture = false;
typeflag.global = false;
typeflag.gradient = false;
addpath('ImFEATbox/GlobalFeatures/Intensity/');
Out = ImFEATbox/GlobalFeatures/Intensity/GradientF(I, typeflag, gradtype);
csvwrite('matlab-out.csv', Out);