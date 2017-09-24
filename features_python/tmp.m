I = csvread('testimg.csv');
addpath('ImFEATbox/GlobalFeatures/Moment/');
Out = ZernikeF(I);
csvwrite('matlab-out.csv', Out);