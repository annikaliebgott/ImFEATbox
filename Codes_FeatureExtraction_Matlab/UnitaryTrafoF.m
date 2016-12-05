function Out = UnitaryTrafoF(I,transformation)
% This function performs various unitary transformations: Walsh-Hadamard,
% Hilbert, Chirp Z, Radon.
%
% Input:     - I: A 2D image
%            - transformation: a struct to determine which transformations
%              should be used for feature extraction.
%              Default: all transformations are carried out
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

if ~exist('transformation','var')
    transformation.WH = true;
    transformation.H = true;
    transformation.R = true;
    transformation.C = true;
end

% convert image
I = double(I);


%% Walsh-Hadamard transform
if transformation.WH
    % Fast Walsh-Hadamard transform
    WH = fwht(I);
    
    %% feature extraction
    % search for strongest coefficients
    WH_max = max(WH(:));
    WH_min = min(WH(:));
    
    % average
    WH_avg = mean2(WH);
    
    % standard deviation
    WH_sd = std2(WH);
    
    % number of nonzero matrix elements, normalized
    WH_nnz = nnz(WH)/numel(WH);
end

%% Hilbert transform
if transformation.H
    % Hilbert transform
    H = hilbert(I);
    
    %% feature extraction
    % search for strongest coefficients
    H_max = max(H(:));
    H_min = min(H(:));
    
    % average
    H_avg = mean2(H);
    
    % standard deviation
    H_sd = std2(H);
    % number of nonzero matrix elements, normalized
    H_nnz  = nnz(H)/numel(H);
end

%% Radon transform
if transformation.R
    % Radon transform
    I = mat2gray(I);
    BW = edge(I);
    theta = 0:179;
    [R,~] = radon(BW,theta);
    
    %% feature extraction
    % search for strongest coefficients
    R_max = max(R(:));
    R_min = min(R(:));
    
    % average
    R_avg = mean2(R);
    
    % standard deviation
    R_sv = std2(R);
    % number of nonzero matrix elements, normalized
    R_nnz = nnz(R)/numel(R);
end

%% Chirp-Z transform
if transformation.C
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
    
    %% feature extraction
    % search for strongest coefficients
    chirp_max = max(chirp(:));
    chirp_min = min(chirp(:));
    
    % chirp matrix = square matrix
    temp1 = eig(chirp);
    chirp_eig = shiftdim(temp1(1:50,1),1);
    chirp_trace = trace(chirp);
    chirp_det = det(chirp);
    chirp_rank = rank(chirp);
    
    % average
    chirp_avg = mean2(chirp);
    
    % standard deviation
    chirp_sv = std2(chirp);
    % number of nonzero matrix elements, normalized
    chirp_nnz = nnz(chirp)/numel(chirp);
end


%% return feature vector

Out = [];
if transformation.WH
    Out = [Out WH_max WH_min WH_avg WH_sd WH_nnz];
end
if transformation.H
    Out = [Out H_max H_min H_avg H_sd H_nnz] ;
end
if transformation.R
    Out = [Out R_max R_min R_avg R_sv R_nnz];
end
if transformation.C
    Out = [Out chirp_max chirp_min chirp_eig chirp_trace...
        chirp_det chirp_rank chirp_avg chirp_sv chirp_nnz];
end



end