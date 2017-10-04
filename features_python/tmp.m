I = csvread('testimg.csv');
addpath('ImFEATbox/LocalFeatures/Point/');
Out = LawF(I);
csvwrite('matlab-out.csv', Out);