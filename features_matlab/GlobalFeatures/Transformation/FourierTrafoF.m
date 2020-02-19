function Out = FourierTrafoF(I,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.global: all features 
%                   + typeflag.transform: all features
%                   + typeflag.moments: only features based on moments
%                   + typeflag.corr: only features based on correlation
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x300) vector containing 300 metrics calculated from
%                   the Fourier transform of an image
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
   typeflag.moments = true;
   typeflag.corr = true;
end    

% change image to a square matrix by means of zero padding
length_x = size(I,1);
length_y = size(I,2);

if (length_x < length_y)
    I = padarray(I, [( length_y - length_x );0],'post');
elseif (length_x > length_y)
    I = padarray(I, [0; ( length_x - length_y )],'post');
elseif (length_x == length_y)
    I = double( I );
end

%reservation for variables
feat_corr = zeros(1,90);
zero_sum = zeros(1,3); zero_sum1 = zeros(1,3); zero_sum2 = zeros(1,3);
zero_sum3 = zeros(1,3); zero_sum4 = zeros(1,3);
power_total = zeros(1,3); power_total1 = zeros(1,3); power_total2 = zeros(1,3);
power_total3 = zeros(1,3); power_total4 = zeros(1,3);
m_2_M = zeros(1,3); m_4_M = zeros(1,3); m_2_M1 = zeros(1,3); m_4_M1 = zeros(1,3);
m_2_M2 = zeros(1,3); m_4_M2 = zeros(1,3); m_2_M3 = zeros(1,3); m_4_M3 = zeros(1,3);
m_2_M4 = zeros(1,3); m_4_M4 = zeros(1,3);
sd_M = zeros(1,3); sd_M1 = zeros(1,3); sd_M2 = zeros(1,3); sd_M3 = zeros(1,3); sd_M4 = zeros(1,3);
sd_F = zeros(1,3); sd_F1 = zeros(1,3); sd_F2 = zeros(1,3); sd_F3 = zeros(1,3); sd_F4 = zeros(1,3);
r_M = zeros(1,3); r_M1 = zeros(1,3); r_M2 = zeros(1,3); r_M3 = zeros(1,3); r_M4 = zeros(1,3);
r_F = zeros(1,3); r_F1 = zeros(1,3); r_F2 = zeros(1,3); r_F3 = zeros(1,3); r_F4 = zeros(1,3);
tr_M = zeros(1,3); tr_M1 = zeros(1,3); tr_M2 = zeros(1,3); tr_M3 = zeros(1,3); tr_M4 = zeros(1,3);
tr_F = zeros(1,3); tr_F1 = zeros(1,3); tr_F2 = zeros(1,3); tr_F3 = zeros(1,3); tr_F4 = zeros(1,3);
max_eig_M = zeros(1,3); max_eig_M1 = zeros(1,3); max_eig_M2 = zeros(1,3);
max_eig_M3 = zeros(1,3); max_eig_M4 = zeros(1,3);
max_eig_F = zeros(1,3); max_eig_F1 = zeros(1,3); max_eig_F2 = zeros(1,3);
max_eig_F3 = zeros(1,3); max_eig_F4 = zeros(1,3);
mean_eig_M = zeros(1,3); mean_eig_M1 = zeros(1,3); mean_eig_M2 = zeros(1,3);
mean_eig_M3 = zeros(1,3); mean_eig_M4 = zeros(1,3);
mean_eig_F = zeros(1,3); mean_eig_F1 = zeros(1,3); mean_eig_F2 = zeros(1,3);
mean_eig_F3 = zeros(1,3); mean_eig_F4 = zeros(1,3);

%% calculate transformation for different orientations

for z=1 : 3
    
    if(z == 2)
        % rotate image by 90 degree
        I = imrotate(I,90,'bilinear','crop');   %'bicubic'
    elseif (z == 3)
        % rotate image by -90 degree
        I = imrotate(I,-180,'bilinear','crop');
    end
    
    % deterimie the Fourier - Transform
    F = fft2(I);
    
    % create translation invariance
    M = abs(mpower(F,2));
    
    % break the image into sectors
    
    % get the center of the image
    half = floor(length(I)/2);
    
    I1 = I(1:half, 1:half);
    I2 = I(1:half,half+1:2*half);
    I3 = I(half+1:2*half, 1:half);
    I4 = I(half+1:2*half, half+1:2*half);
    
    % Fourier Transform
    F1 = fft2(I1);
    F2 = fft2(I2);
    F3 = fft2(I3);
    F4 = fft2(I4);
    
    % create translation invariance
    M1 = abs((F1.^2));
    M2 = abs((F2.^2));
    M3 = abs(F3.^2);
    M4 = abs((F4.^2));
    
    %% extract features
    
    % determine the moments
    if (typeflag.moments || typeflag.global || typeflag.transform)
        m_2_M(z) = mean2(moment(M,2));
        m_4_M(z) = mean2(moment(M,4));
        m_2_M1(z) = mean2(moment(M1,2));
        m_4_M1(z) = mean2(moment(M1,4));
        m_2_M2(z) = mean2(moment(M2,2));
        m_4_M2(z) = mean2(moment(M2,4));
        m_2_M3(z) = mean2(moment(M3,2));
        m_4_M3(z) = mean2(moment(M3,4));
        m_2_M4(z) = mean2(moment(M4,2));
        m_4_M4(z) = mean2(moment(M4,4));
    end
    
    if (typeflag.global || typeflag.transform)
        % standard derivation
        sd_M(z) = std2(M);
        sd_M1(z) = std2(M1);
        sd_M2(z) = std2(M2);
        sd_M3(z) = std2(M3);
        sd_M4(z) = std2(M4);
        
        sd_F(z) = std2(F);
        sd_F1(z) = std2(F1);
        sd_F2(z) = std2(F2);
        sd_F3(z) = std2(F3);
        sd_F4(z) = std2(F4);
        
        % rank
        r_M(z) = rank(M);
        r_M1(z) = rank(M1);
        r_M2(z) = rank(M2);
        r_M3(z) = rank(M3);
        r_M4(z) = rank(M4);
        
        r_F(z) = rank(F);
        r_F1(z) = rank(F1);
        r_F2(z) = rank(F2);
        r_F3(z) = rank(F3);
        r_F4(z) = rank(F4);
        
        % trace
        tr_M(z) = trace(M);
        tr_M1(z) = trace(M1);
        tr_M2(z) = trace(M2);
        tr_M3(z) = trace(M3);
        tr_M4(z) = trace(M4);
        
        tr_F(z) = trace(F);
        tr_F1(z) = trace(F1);
        tr_F2(z) = trace(F2);
        tr_F3(z) = trace(F3);
        tr_F4(z) = trace(F4);
        
        % largest/smallest eingenvalue
        eig_M = eig(M);
        eig_M1 = eig(M1);
        eig_M2 = eig(M2);
        eig_M3 = eig(M3);
        eig_M4 = eig(M4);
        
        eig_F = eig(F);
        eig_F1 = eig(F1);
        eig_F2 = eig(F2);
        eig_F3 = eig(F3);
        eig_F4 = eig(F4);
        
        max_eig_M(z) = max(eig_M);
        max_eig_M1(z) = max(eig_M1);
        max_eig_M2(z) = max(eig_M2);
        max_eig_M3(z) = max(eig_M3);
        max_eig_M4(z) = max(eig_M4);
        
        max_eig_F(z) = max(eig_F);
        max_eig_F1(z) = max(eig_F1);
        max_eig_F2(z) = max(eig_F2);
        max_eig_F3(z) = max(eig_F3);
        max_eig_F4(z) = max(eig_F4);
        
        % mean2
        mean_eig_M(z) = mean2(eig_M);
        mean_eig_M1(z) = mean2(eig_M1);
        mean_eig_M2(z) = mean2(eig_M2);
        mean_eig_M3(z) = mean2(eig_M3);
        mean_eig_M4(z) = mean2(eig_M4);
        
        mean_eig_F(z) = mean2(eig_F);
        mean_eig_F1(z) = mean2(eig_F1);
        mean_eig_F2(z) = mean2(eig_F2);
        mean_eig_F3(z) = mean2(eig_F3);
        mean_eig_F4(z) = mean2(eig_F4);
    end
    
    % correlation
    if (typeflag.corr || typeflag.transform || typeflag.global)
        corr12 = corr(M1,M2);
        corr13 = corr(M1,M3);
        corr14 = corr(M1,M4);
        corr32 = corr(M3,M2);
        corr42 = corr(M4,M2);
        corr34 = corr(M4,M3);
        feat_corr((z-1)*30+1:(z-1)*30+30) = [std2(corr12) mean2(corr12)...
            max(corr12(:)) min(corr12(:)) nnz(corr12)...
            std2(corr13) mean2(corr13) max(corr13(:)) min(corr13(:)) nnz(corr13)...
            std2(corr14) mean2(corr14) max(corr14(:)) min(corr14(:)) nnz(corr14)...
            std2(corr32) mean2(corr32) max(corr32(:)) min(corr32(:)) nnz(corr32)...
            std2(corr42) mean2(corr42) max(corr42(:)) min(corr42(:)) nnz(corr42)...
            std2(corr34) mean2(corr34) max(corr34(:)) min(corr34(:)) nnz(corr34)];
    end
    
    if (typeflag.global || typeflag.transform)
        % count zero points
        zero_sum(z) = nnz(F);
        zero_sum1(z) = nnz(F1);
        zero_sum2(z) = nnz(F2);
        zero_sum3(z) = nnz(F3);
        zero_sum4(z) = nnz(F4);
        
        
        % compute the DFT power
        
        % Transform length
        n = pow2(nextpow2(length(I)));
        
        % Power of the DFT
        power = F.*conj(F);
        power1 = F1.*conj(F1);
        power2 = F2.*conj(F2);
        power3 = F3.*conj(F3);
        power4 = F4.*conj(F4);
        power_total(z) = sum( sum(power) ) /  (n*n);
        power_total1(z) = sum( sum(power1) ) / ( n*n );
        power_total2(z) = sum( sum(power2) ) / ( n*n );
        power_total3(z) = sum( sum(power3) ) / ( n*n );
        power_total4(z) = sum( sum(power4) ) / ( n*n );
    end
    
    
end

%% return feature vector
if (typeflag.global || typeflag.transform)
    Out = [zero_sum zero_sum1 zero_sum2 zero_sum3 zero_sum4...
        power_total power_total1 power_total2 power_total3 power_total4...
        m_2_M m_4_M m_2_M1 m_4_M1 m_2_M2 m_4_M2 m_2_M3 m_4_M3 m_2_M4 m_4_M4...
        sd_M sd_M1 sd_M2 sd_M3 sd_M4...
        sd_F sd_F1 sd_F2 sd_F3 sd_F4...
        r_M r_M1 r_M2 r_M3 r_M4...
        r_F r_F1 r_F2 r_F3 r_F4...
        tr_M tr_M1 tr_M2 tr_M3 tr_M4...
        tr_F tr_F1 tr_F2 tr_F3 tr_F4...
        max_eig_M max_eig_M1 max_eig_M2 max_eig_M3 max_eig_M4...
        max_eig_F max_eig_F1 max_eig_F2 max_eig_F3 max_eig_F4...
        mean_eig_M mean_eig_M1 mean_eig_M2 mean_eig_M3 mean_eig_M4...
        mean_eig_F mean_eig_F1 mean_eig_F2 mean_eig_F3 mean_eig_F4...
        feat_corr];
elseif (~typeflag.moments)
    Out = feat_corr;
elseif (~typeflag.corr)
    Out = [m_2_M m_4_M m_2_M1 m_4_M1 m_2_M2 m_4_M2 m_2_M3 m_4_M3 m_2_M4 m_4_M4];
else
    Out = [m_2_M m_4_M m_2_M1 m_4_M1 m_2_M2 m_4_M2 m_2_M3 m_4_M3 m_2_M4 m_4_M4...
        feat_corr];
end

end
