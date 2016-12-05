function [Out] = GLCMF(Image,InputParameters, typeflag)
% Input:     - I: A 2D image
%            - InputParameters : A structure containing parameters:
%                   + DisplacementVector: A (nx2) vector composed of offset 
%                   and orientation (default = [0 1])
%                   + NumLevels: An integer number specifying the number of 
%                   gray level intensities (default = 8)
%                   + GrayLimits: A (1x2) vector containing the min/max 
%                   gray levels which are needed to sort gray levels of I
%                   into number of gray levels specified by NumLevels
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.global: all features 
%                   + typeflag.form: all features
%                   + typeflag.texture: all features
%                   + typeflag.corr: only features based on correlation
%                   + typeflag.entropy: only features based on entropy
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1xn*21) vector containing n*21 metrics calculated 
%              from the gray level co-occurence matrix
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
% Implementation based on:  R. Haralick; K. Shanmugam; I. Dinstein (1973): 
%                           "Textural Features for Image Classification". 
%                           IEEE Transactions on Systems, Man, and 
%                           Cybernetics. SMC-3 (6): 610–621.


% Check for color image and convert to the grayscale
if(numel(size(Image))==3)       %if image is 3D
    if(size(Image,3)==3)
        Image = rgb2gray(Image);
    end
end


% Default Parameters (Offset value = 1 pixel)
if ~exist('InputParameters','var') || ~isfield(InputParameters,'DisplacementVector') 
    InputParameters.DisplacementVector = [0 1];
end

if ~exist('InputParameters','var') || ~isfield(InputParameters,'NumLevels')
    InputParameters.NumLevels = 8;
end

if ~exist('InputParameters','var') || ~isfield(InputParameters,'GrayLimits')
    InputParameters.GrayLimits = [min(Image(:)) max(Image(:))];
end

if ~exist('typeflag', 'var')
   typeflag.global = true; 
   typeflag.texture = true;
   typeflag.form = true;
   typeflag.corr = true;
   typeflag.entropy = true;
end    

if typeflag.global || typeflag.texture
    typeflag.corr = true;
    typeflag.entropy = true;
end    

GLCM_Matrices = graycomatrix(Image, 'offset', InputParameters.DisplacementVector,...
    'NumLevels',InputParameters.NumLevels,...
    'GrayLimits', InputParameters.GrayLimits,...
    'Symmetric', false);

if typeflag.texture || typeflag.global
    % extract all features of GLCM
    Out = zeros(size(GLCM_Matrices,3),21);
elseif typeflag.corr
    if typeflag.entropy
        % ectract correlation based and entropy based features
        Out = zeros(size(GLCM_Matrices,3),7);
    else
        % extract only correlation based features
        Out = zeros(size(GLCM_Matrices,3),4);
    end
else
    % extract only entropy based features
    Out = zeros(size(GLCM_Matrices,3),5);
end


for n=1:1:size(GLCM_Matrices,3)
    
    G = GLCM_Matrices(:,:,n);
    SumOfGLCM = sum(sum(G));
    GNormalized = G ./ SumOfGLCM;
    MeanOfGLCM = mean2(GNormalized);
    s1 = size(G,1);
    s2 = size(G,2);
    
    % Initialization
    p_x   = zeros(s1,1);
    p_y   = zeros(s2,1);
    u_x   = 0;
    u_y   = 0;
    HXY1  = 0;
    HXY2  = 0;
    sigma_x = 0;
    sigma_y = 0;
    p_xplusy = zeros((s1*2 - 1),1);
    p_xminusy = zeros(s1,1);
    CO      = 0;
    CORR    = 0;
    IDM     = 0;
    SE      = 0;
    H       = 0;
    DE      = 0;
    IMC1    = 0;
    IMC2    = 0;
    ACORR   = 0;
    DIS     = 0;
    CLS     = 0;
    CLP     = 0;
    INV     = 0;
    INVN     = 0;
    IDMN     = 0;
    
    
    m = 1:s1; %replace one for-loop by using a vector of length s1
    o = 1:s2; %replace one for-loop by using a vector of length s2
    
    % Angular Second Moment or Energy (ASM)
    ASM = sum(GNormalized(:).^2);
    
    % Sum of squared variance (SSV)
    SSV = sum(((m - MeanOfGLCM).^2)*GNormalized);
    
    % Entropy (H)
    if typeflag.entropy
        H = - sum(GNormalized(:).*log(GNormalized(:)));
    end
    
    for i = 1:s1
        
        u_x = u_x + sum(i*GNormalized(i,:));
        u_y = u_y + sum(o.*GNormalized(i,:));
        
        p_x(i) = sum(GNormalized(i,:));
        p_y(i) = sum(GNormalized(:,i));
        
        k = i+(1:s2);
        l = abs(i-(1:s2));
        p_xplusy(k-1) = p_xplusy(k-1) + GNormalized(i,(1:s2)).';
        
        % since there are multiple equal values in index vector l (i.e. some
        % indices are in there twice) and matlab can't automatically assign
        % the sum of multiple values to the same index, the calculation has
        % to be divided into two parts (on from first index down until l==1,
        % then up from l==0 to l(end)
        p_xminusy(l(1:i-1)+1) = p_xminusy(l(1:i-1)+1) + GNormalized(i,(1:i-1)).';
        p_xminusy(l(i:end)+1) = p_xminusy(l(i:end)+1) + GNormalized(i,(i:s2)).';
        
        
        % Contrast (CO)
        CO = CO + sum((i - o).^2.*GNormalized(i,:));
        
        % Inverse difference moment or homogenity (IDM)
        IDM = IDM + sum(GNormalized(i,o)./( 1 + (i - o).^2));
        
        % Dissimilarity (DIS)
        DIS = DIS + sum(abs(i - o).*GNormalized(i,:));
        
        % Autocorrelation (ACORR)
        if typeflag.corr
            ACORR = ACORR + sum(i*o.*GNormalized(i,o));
        end
        
        % Inverse difference (INV)
        INV = INV + sum(GNormalized(i,o)./( 1 + abs(i - o)));
        
        % Inverse difference normalized (INVN)
        INVN = INVN + sum(GNormalized(i,o)./( 1 + (abs(i-o)./(s1)) ));
        
        % Inverse difference moment normalized (IDMN)
        IDMN = IDMN + sum(GNormalized(i,o)./( 1 + ((i-o)./(s1)).^2));
        
    end
    
    for i = 1:s1
        
        % Cluster shade (CLS)
        CLS = CLS + ((i + o - u_x - u_y).^3)*GNormalized(i,o)';
        
        % Cluster prominence (CLP)
        CLP = CLP + ((i + o - u_x - u_y).^4)*GNormalized(i,o)';
        
        sigma_x  = sigma_x  + sum((((i) - u_x)^2).*GNormalized(i,o));
        sigma_y  = sigma_y  + sum((((o) - u_y).^2).*GNormalized(i,o));
        
        
    end
    
    
    % Correlation (CORR)
    if typeflag.corr
        CORR = (ACORR - (u_x*u_y)) / (sqrt(sigma_x)*sqrt(sigma_y));
    end
    
    % Summed average (SA)
    SA =(2:2*s1)*p_xplusy;
    
    % Summed entropy (SE)
    if typeflag.entropy
        SE = -sum(p_xplusy.*log(p_xplusy));
    end
    
    % Summed Variance (SV)
    SV = (((2:2*s1) - SE).^2)*p_xplusy;
    
    % Difference varience (DV)
    DV = ((0:s1-1).^2)*p_xminusy;
    
    % Difference entropy (DE)
    if typeflag.entropy
        DE = -sum(p_xminusy.*log(p_xminusy));
        
        HXY = H;
        
        for i = 1:s1
            HXY1 = HXY1 - GNormalized(i,:)*log(p_x(i).*p_y);
            HXY2 = HXY2 - sum(p_x(i).*p_y.*log(p_x(i).*p_y));
        end
        
        Hx = - sum(p_x.*log(p_x));
        Hy = - sum(p_y.*log(p_y));
        
        % Information measure of correlation 1 (IMC1)
        % Information measure of correlation 2 (IMC2)
        IMC1 = ( HXY - HXY1) / ( max([Hx,Hy]) );
        IMC2 = ( 1 - exp(-2*(HXY2 - HXY)))^0.5;
    end
    
    % Maximum probability (MAXP)
    MAXP = max(GNormalized(:));
    
    if typeflag.texture || typeflag.global
        Out(n,:) = [ACORR , CO, CORR, CLP, CLS ,  DIS, ASM, H, IDM, MAXP,...
            SSV, SA, SV, SE, DV, DE, IMC1, IMC2, INV,INVN, IDMN];
    elseif typeflag.corr
        if typeflag.entropy
            Out(n,:) = [ACORR, CORR, H, SE, DE, IMC1, IMC2];
        else
            Out(n,:) = [ACORR, CORR, IMC1, IMC2];
        end
    else
        Out(n,:) = [H, SE, DE, IMC1, IMC2];
    end
end

Out = reshape(Out, [1,size(Out,1)*size(Out,2)]);

end

