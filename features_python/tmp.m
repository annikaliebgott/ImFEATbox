I = csvread('testimg.csv');
addpath('ImFEATbox/GlobalFeatures/Geometrical/');
Out = RunLengthF(I);
csvwrite('matlab-out.csv', Out);