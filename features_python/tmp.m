I = csvread('testimg.csv');
plotflag = false;
addpath('ImFEATbox/LocalFeatures/');
Out = HarrisF(I, plotflag);
csvwrite('matlab-out.csv', Out);