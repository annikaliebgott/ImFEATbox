I = csvread('testimg.csv');
addpath('ImFEATbox/GlobalFeatures/Moment/');
Out = HuF(I);
csvwrite('matlab-out.csv', Out);