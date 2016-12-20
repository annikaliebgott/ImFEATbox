function Out = DistanceTrafoF(I,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.global: all features 
%                   + typeflag.transform: all features 
%                   + typeflag.corr: only features based on correlation
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x56) vector containing 56 metrics calculated from
%              the distance transform
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic 
% and Interventional Radiology, University Hospital of Tuebingen, Germany 
% and the Institute of Signal Processing and System Theory University of 
% Stuttgart, Germany. Last modified: November 2016
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
% ************************************************************************

if ~exist('typeflag','var')
   typeflag.global = true; 
   typeflag.transform = true;
   typeflag.corr = true;
end    

% converte image
% im2bw can't process complex input values
I = double(real(I)); 

% graythresh(): Global image threshold using Otsu's method
BW = im2bw(I, graythresh(I));

%% transfom image
[D1, IDX1] = bwdist(BW, 'chessboard');
[D2, IDX2] = bwdist(BW, 'cityblock');
[D3, IDX3] = bwdist(BW, 'euclidean');
[D4, IDX4] = bwdist(BW, 'quasi-euclidean');


%% feature extraction

% 2D correlation coefficient
r12 = corr2(IDX1,IDX2);
r13 = corr2(IDX1,IDX3);
r14 = corr2(IDX1,IDX4);
r23 = corr2(IDX2,IDX3);
r24 = corr2(IDX2,IDX4);
r34 = corr2(IDX3,IDX4);

rD12 = corr2(D1,D2);
rD13 = corr2(D1,D3);
rD14 = corr2(D1,D4);
rD23 = corr2(D2,D3);
rD24 = corr2(D2,D4);
rD34 = corr2(D3,D4);


if (typeflag.global || typeflag.transform)
    % mean of matrix elements
    B1 = mean2(IDX1);
    B2 = mean2(IDX2);
    B3 = mean2(IDX3);
    B4 = mean2(IDX4);
    
    B11 = mean2(D1);
    B22 = mean2(D2);
    B33 = mean2(D3);
    B44 = mean2(D4);
    
    % standard deviation of matrix elements
    S1 = std2(IDX1);
    S2 = std2(IDX2);
    S3 = std2(IDX3);
    S4 = std2(IDX4);
    
    S11 = std2(D1);
    S22 = std2(D2);
    S33 = std2(D3);
    S44 = std2(D4);
    
    s1 = std(std(double(IDX1)));
    s2 = std(std(double(IDX2)));
    s3 = std(std(double(IDX3)));
    s4 = std(std(double(IDX4)));
    
    s11 = std(std(double(D1)));
    s22 = std(std(double(D2)));
    s33 = std(std(double(D3)));
    s44 = std(std(double(D4)));
    
    % number of non zero elements
    nn_D11 = nnz(IDX1);
    nn_D22 = nnz(IDX2);
    nn_D33 = nnz(IDX3);
    nn_D44 = nnz(IDX4);
    
    nn_D1 = nnz(D1);
    nn_D2 = nnz(D2);
    nn_D3 = nnz(D3);
    nn_D4 = nnz(D4);
    
    % determine max/min values of IDX
    max_IDX1 = max(IDX1(:));
    max_IDX2 = max(IDX2(:));
    max_IDX3 = max(IDX3(:));
    max_IDX4 = max(IDX4(:));
    min_IDX1 = min(IDX1(:));
    min_IDX2 = min(IDX2(:));
    min_IDX3 = min(IDX3(:));
    min_IDX4 = min(IDX4(:));
    
    % determine max value of D
    max_D1 = max(D1(:));
    max_D2 = max(D2(:));
    max_D3 = max(D3(:));
    max_D4 = max(D4(:));
    
end

%% return feature vector

if ~(typeflag.global || typeflag.transform)
    Out = [r12 r13 r14 r23 r24 r34 rD12 rD13 rD14 rD23 rD24 rD34];
else
    Out = [B1 B2 B3 B4 B11 B22 B33 B44...
        S1 S2 S3 S4 S11 S22 S33 S44...
        s1 s2 s3 s4 s11 s22 s33 s44...
        r12 r13 r14 r23 r24 r34 rD12 rD13 rD14 rD23 rD24 rD34...
        nn_D1 nn_D2 nn_D3 nn_D4 nn_D11 nn_D22 nn_D33 nn_D44...
        max_IDX1 max_IDX2 max_IDX3 max_IDX4 min_IDX1 min_IDX2 min_IDX3 min_IDX4...
        max_D1 max_D2 max_D3 max_D4];
end

end
