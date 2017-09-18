I = csvread('testimg.csv');
gradtype.first = true;
gradtype.second = true;
typeflag.texture = true;
typeflag.global = true;
typeflag.gradient = true;
addpath('ImFEATbox/GlobalFeatures/Intensity/');
Out = GradientF(I, typeflag, gradtype);
csvwrite('matlab-out.csv', Out);