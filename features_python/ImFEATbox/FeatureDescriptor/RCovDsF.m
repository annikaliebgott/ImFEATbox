function Out = RCovDsF(I,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.corr: all features 
%                   + typeflag.moments: only features based on moments
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x15) vector containing 15 metrics calculated based
%              on Region Covariance Descriptors
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
%
% Implementation based on:  M. Faraki, M. T. Harandi and F. Porikli, 
%                           "Approximate infinite-dimensional Region 
%                           Covariance Descriptors for image 
%                           classification," 2015 IEEE International 
%                           Conference on Acoustics, Speech and Signal
%                           Processing (ICASSP), South Brisbane, QLD, 2015, 
%                           pp. 1364-1368.

if ~exist('typeflag','var')
   typeflag.corr = true;
   typeflag.moments = true;
end    

% convert image to double
I = double(I);

% dimensions of I
n = size(I,2);
d = size(I,1);


%% Random Fourier Features
% arbitrary kernel is not shift invariant
% Sample mean of the observations

% define matrix J
J = sqrt(mpower(n,-3))*(n*eye(n)-ones(1,n)*ones(n,1));

% i.i.d samples drawn from N(0,(sigma^-2)*I)
omega_i = double(samplesN(d,n));

% uniform samples drawn from [0 2pi]
q = 2*pi/size(I,2);
b_i = double(0:q:2*pi);

% mapping
z_F = zeros(1,size(I,2));
for i=1 : size(I,2)
    z_F(i) = sqrt(2/d)*cos(transpose(omega_i(:,i))*I(:,i) + b_i(i));
end

RCovD_fourier = z_F*J*transpose(J)*transpose(z_F);

% extract 2nd and 4th moment as features
z_F_moment2 = moment(z_F,2);
z_F_moment4 = moment(z_F,4);


%% Nystroem Method
% data dependent estimation of Reducing Kernel Hilbert Space (RKHS) by
% using the Nystroem Mehtod

% collection of training examples
rng('default')
D = I(:,randperm(size(I,2),90));

% compute the Kernel matrix K, ker(A) = u with A*u = 0
% ZR = null(D,'r');
sig = std2(D);
if sig ~= 0
    n2 = size(D,1);
    dia = diag(D*D'/sig^2);
    K = exp(D*D'/(sig^2)-ones(n2,1)*dia'/2-dia*ones(1,n2)/2);
    
    % diagonal matrix of top D eigenvalues of K
    [V, Sigma] = eig(K);
    
    
    % associated eigenvalues
    z_N = zeros(size(D,1),size(I,2));
    for j=1 : size(I,2)
        k_x = zeros(1,size(D,2));
        for k=1 : size(D,2)
            k_x(k) = mean((2/d)*cos(transpose(omega_i(:,k))*D(:,k)+b_i(k))...
                *cos(transpose(omega_i(:,j))*I(:,j)+b_i(j)));
        end
        % fill k_x with zeros to fit the dimensions
        k_x = padarray(k_x, [0 (size(D,1)-length(k_x))],'post');
        z_N(:,j) = (sqrtm(Sigma))\(V*transpose(k_x));
    end
    
    
    % compute RCovD_nystr
    RCovD_nystr = z_N*J*transpose(J)*transpose(z_N);
else
    RCovD_nystr = 0;
    z_N = 0;
end

% extract features from moments
RCovD_nystr_moment2 = moment(z_N,2);
RCovD_nystr_moment4 = moment(z_N,4);
RCovD_nystr_moment2_mean = mean(RCovD_nystr_moment2);
RCovD_nystr_moment4_mean = mean(RCovD_nystr_moment4);
RCovD_nystr_moment2_std = std(RCovD_nystr_moment2);
RCovD_nystr_moment4_std = std(RCovD_nystr_moment4);
RCovD_nystr_moment2_min = min(RCovD_nystr_moment2);
RCovD_nystr_moment4_min = min(RCovD_nystr_moment4);
RCovD_nystr_moment2_max = max(RCovD_nystr_moment2);
RCovD_nystr_moment4_max = max(RCovD_nystr_moment4);

if typeflag.corr
    RCovD_nystr_mean = mean2(RCovD_nystr);
    RCovD_nystr_std = std2(RCovD_nystr);
    RCovD_nystr_max = max(RCovD_nystr(:));
    RCovD_nystr_min = min(RCovD_nystr(:));
end

%% return feature vector

if typeflag.corr
    Out = [RCovD_fourier z_F_moment2 z_F_moment4...
        RCovD_nystr_moment2_mean RCovD_nystr_moment4_mean...
        RCovD_nystr_moment2_std RCovD_nystr_moment4_std...
        RCovD_nystr_moment2_min RCovD_nystr_moment4_min...
        RCovD_nystr_moment2_max RCovD_nystr_moment4_max...
        RCovD_nystr_mean RCovD_nystr_std RCovD_nystr_max RCovD_nystr_min];
else
    Out =  [z_F_moment2 z_F_moment4...
        RCovD_nystr_moment2_mean RCovD_nystr_moment4_mean...
        RCovD_nystr_moment2_std RCovD_nystr_moment4_std...
        RCovD_nystr_moment2_min RCovD_nystr_moment4_min...
        RCovD_nystr_moment2_max RCovD_nystr_moment4_max];
end

end

function z = samplesN(d,n)
%return a matrix which has normally distributed enties

npdf = load('normaly_distributed_matrix.mat');
npdf = cell2mat(struct2cell(npdf));
z = npdf(1:d,1:n);
end



