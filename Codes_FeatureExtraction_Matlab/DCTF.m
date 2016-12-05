function Out = DCTF(I,typeflag)
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
% Output:    - Out: A (1x2901) vector containing 2901 metrics calculated
%                   from the discrete cosine transform
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

I = double(I);

% reserve space for the feature vector
% if typeflag.transform || typeflag.global == true: extract all features
% if only typeflag.corr == true: extract especially correlation based 
% features
if typeflag.transform || typeflag.global
    f = zeros(3,967);
else
    f = zeros(3,84);
end

%% calculate DCT for different orientations: horizontal,vertical, diagonal

for z = 1:3
    % horizontal, vertical and diagonal by counting the number of occurrence
    % of each selected coefficient with respect to its position
    % equal importance of horizontal, vertical and diagonal
    
    if (z == 2)
        % rotate image by 90 degree
        I = imrotate(I,90,'bilinear','crop');   %'bicubic'
    elseif (z == 3)
        % rotate image by -90 degree
        I = imrotate(I,-180,'bilinear','crop');
    end
    
    % index counter
    idx = 0;
    
    % 2D discrete cosine transform, later on used as reference matrix
    B = dct(I);
    
    % perform SVD
    [U,S,V] = svd(B);
    eUB = eig(U);
    eVB = eig(V);
    
    % extract features
    if typeflag.transform || typeflag.global
        f(z,idx+1:idx+100) = [shiftdim(eUB(1:50), 1) shiftdim(eVB(1:50), 1)];
        idx = idx+100;

        coefB = [B(1,1), B(1,40), B(20,60), B(80,100) B(100, 150)];
        f(z,idx+1:idx+length(coefB)) = coefB;
        idx = idx + length(coefB);

        f(z,idx+1:idx+4) = [det(U) trace(U) det(V) trace(V)];
        idx = idx + 4;
    end
    
    %% calculate transform for different block decomposition sizes
    for i = [2, 4, 8, 16, 32, 64]
        
        % extract spatio-temporal features
        % decompose blocks of images into independent tiles
        fun = @(block_struct) ...
            std2(block_struct.data) * ones(size(block_struct.data));
        I2 = blockproc(I,[i i],fun);
        
        
        % 2D discrete cosine transform
        % I2 has to be grayscale
        B2 = dct(I2);
        
        % output of some coefficients
        % 3D zigzag transversal to select 25% of the coefficents
        if typeflag.transform || typeflag.global        
            coefB2 = [B(1,1), B(1,40), B(20,60), B(80,100) B(100, 150)];
            f(z,idx+1:idx+length(coefB2)) = coefB2;
            idx = idx + length(coefB2);

            % calculate important features of the 2D DCT
            temp = [std2(B2) std(std(B2)) mean2(B2) rank(B2) max(B2(:))...
                min(B2(:)) nnz(B2)];
            f(z,idx+1:idx+length(temp)) = temp;
            idx = idx + length(temp);
        end
        
        % correlation between DCT of I and decomposed I
        Corr = corrcoef(B, B2);
        f(z,idx+1) = Corr(1,2);
        idx = idx +1;
        
        % singular value deomposition, U & V are square matrices
        [U2,S2,V2] = svd(B2);
        eUB2 = eig(U2);
        eVB2 = eig(V2);
        
        if typeflag.transform || typeflag.global
            temp = [shiftdim(eUB2(1:50),1) shiftdim(eVB2(1:50), 1)...
                std2(U2) std2(S2) std2(V2) mean2(U2) mean2(S2) mean2(V2)...
                max(U2(:)) max(S2(:)) max(V2(:)) min(U2(:)) min(S2(:)) min(V2(:))...
                nnz(B2) det(U2) det(V2) trace(U2) trace(V2)];
            f(z,idx+1:idx+length(temp)) = temp;
            idx = idx + length(temp);
        end 
        
        % correlation between SVD matrices
        corrS = corrcoef(S,S2);
        corrU = corrcoef(U,U2);
        corrV = corrcoef(V,V2);
        
        % correlation between eigenvalues of SVD
        correV = corrcoef(eVB,eVB2);
        correU = corrcoef(eUB,eUB2);
        xcorreV = xcorr(eVB,eVB2);
        xcorreU = xcorr(eUB,eUB2);
        temp = [corrS(1,2) corrU(1,2) corrV(1,2) correV(1,2) correU(1,2)...
            mean(xcorreV) std(xcorreV) min(xcorreV) max(xcorreV)...
            mean(xcorreU) std(xcorreU) min(xcorreU) max(xcorreU)];
        
        f(z,idx+1:idx+length(temp)) = temp;
        idx = idx + length(temp);
    end
end

%% return feature vector
f = shiftdim(f,1);
Out = f(:).';

end