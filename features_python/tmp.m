I = csvread('testimg.csv');
typeflag.local = true;
typeflag.moments = true;
typeflag.corr = true;
typeflag.texture = true;
addpath('ImFEATbox/LocalFeatures/Line/');
Out = ImFEATbox/LocalFeatures/Line/LineProfileF(I, plotflag, typeflag);
csvwrite('matlab-out.csv', Out);