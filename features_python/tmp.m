I = csvread('testimg.csv');
plotflag = false;
typeflag.moments = true;
typeflag.local = true;
typeflag.corr = true;
typeflag.texture = true;
addpath('ImFEATbox/LocalFeatures/Line/');
Out = LineProfileF(I, plotflag, typeflag);
csvwrite('matlab-out.csv', Out);