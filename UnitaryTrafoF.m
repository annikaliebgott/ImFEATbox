function Out = UnitaryTrafoF(I)
% This function performs various unitary transformations: Walsh-Hadamard,
% Hilbert, Chirp Z, Radon.
%
% Input:     - I: A 2D image
%
%
% Output:    - Out: A (1x73) vector containing 73 metrics calculated from
%              various unitary transformations of the image
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


%% Transformation

% convert image
I = double(I);

% Fast Walsh-Hadamard transform
WH = fwht(I);

% Hilbert transform
H = hilbert(I);

% form matrix to sqaure matirx
B=I;
s = size(B);
if (s(1,1) < s(1,2))
    B = padarray(B, [( s(1,2) - s(1,1) );0],'post');
elseif (s(1,1) > s(1,2))
    B = padarray(B, [0; ( s(1,1) - s(1,2) )],'post');
elseif (s(1,1) == s(1,2))
    B = double( B );
end


% Chirp Z-transform
chirp = czt(B);

% Radon transform
I = mat2gray(I);
BW = edge(I);
theta = 0:179;
[R,~] = radon(BW,theta);


%% feature extraction

% search for strongest coefficients
WH_max = max(WH(:));
WH_min = min(WH(:));
H_max = max(H(:));
H_min = min(H(:));
chirp_max = max(chirp(:));
chirp_min = min(chirp(:));
R_max = max(R(:));
R_min = min(R(:));

% chirp matrix = square matrix
temp1 = eig(chirp);
chirp_eig = shiftdim(temp1(1:50,1),1);
chirp_trace = trace(chirp);
chirp_det = det(chirp);
chirp_rank = rank(chirp);

% average
WH_avg = mean2(WH);
H_avg = mean2(WH);
chirp_avg = mean2(chirp);
R_avg = mean2(R);

% standard deviation
WH_sd = std2(WH);
H_sd = std2(H);
chirp_sv = std2(chirp);
R_sv = std2(R);

% number of nonzero matrix elements, normalized
WH_nnz = nnz(WH)/numel(WH);
H_nnz  = nnz(H)/numel(H);
R_nnz = nnz(R)/numel(R);
chirp_nnz = nnz(chirp)/numel(chirp);




%% return feature vector
Out = [WH_max WH_min WH_avg WH_sd WH_nnz...
    H_max H_min H_avg H_sd H_nnz...
    chirp_max chirp_min chirp_eig chirp_trace...
    chirp_det chirp_rank chirp_avg chirp_sv chirp_nnz...
    R_max R_min R_avg R_sv R_nnz
    ];

end