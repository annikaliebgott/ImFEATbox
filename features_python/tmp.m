I = csvread('testimg.csv');
typeflag.form = true;
typeflag.global = true;
addpath('ImFEATbox/GlobalFeatures/Geometrical/');
Out = FormFactorF(I, typeflag);
csvwrite('matlab-out.csv', Out);