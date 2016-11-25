%% Main script to carry out the feature extraction
%
% This script calls all feature extraction allgorithms provided by
% ImFEATbox. 
%
% Output:   feat_vec: an [NxM] array of M features extracted from N images.
%
% The basic steps to successfully carry out the feature extraction:
% 1.)   Import the images of which you want to extract the features. 
%       Note: this script expects images to be cell arrays
% 2.)   Set typeflag to define which features you wish to extract.
%       For more information about typeflag: see README.txt
% 3.)   Set parameters needed for some feature extraction algorithms.
%       Note: default values are tuned for use with automated MR image
%       quality assessment, you might need to change them according to your
%       application
% 4.)   Choose wheather or not you wish to use additional tools
%       (preprocessing, visualization)
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


%% import images to extract the features from
% Note: this script expects images to be saved as cell arrays
data = load('images.mat');
images = data.images;

% Determine total number of 2D images to extract features from
N_slices_total = 0;
N_images = length(images);
for i = 1:length(images)
    N_slices_total = N_slices_total + size(images{i},3);
end

%% set typeflag to choose which types of features are being extracted
typeflag = struct;

% If you wish to extract all features provided by ImFEATbox, set
% typeflag.all = true
typeflag.all = false;

% If you don't wish to extract all features, you can chose the desired
% feature categories here
% Note: if typeflag.all = true, all other flags are automatically also set 
% to true. All changes set manually will be overwritten!

typeflag.global = true;
typeflag.local = false;
typeflag.corr = false;
typeflag.gradient = false;
typeflag.moments = false;
typeflag.texture = false;
typeflag.form = false;
typeflag.entropy = false;
typeflag.transform = false;

if typeflag.all
    typeflag.global = true;
    typeflag.local = true;
    typeflag.corr = true;
    typeflag.gradient = true;
    typeflag.moments = true;
    typeflag.texture = true;
    typeflag.form = true;
    typeflag.entropy = true;
    typeflag.transform = true;
end
 

%% Set parameters
% set GLCM parameters for GLCM feature extraction
values = [1;2;4;8;16;32;64;128];
offset0 = [zeros(size(values,1),1) values];
offset1 = [-values values];
offset2 = [-values zeros(size(values,1),1)];
offset3 = [-values -values];
offset = [offset0 ; offset1 ; offset2 ; offset3];
GLCMParameters.DisplacementVector = offset;
GLCMParameters.NumLevels = 255;

% set filter parameters for LBP feature extraction
RadialLBP = [8 1; 16 2;24 3;32 4];

% set desired number of blobs for Laplacian of Gaussian feature extraction
N_blobs = 120;

% set threshold for Quadtree decomposition feature extraction
threshold_quadtree = 0.27;

% set parameters for Hough transform features
houghtype = 'both';
arc_min = pi/2;

%% Choose optional steps
% set which preprocessing steps are desired
doSegmentation = true;
doGrayscaling = true;

% set whether visualizations included in some feature extraction algorithms
% should be plotted or not (note: increases computation time)
plotflag = false;


%% Preallocate feature arrays for speed 
% implemented only for settings where a high number of features are being 
% extracted: all features, all global features, all local features
if typeflag.all
    feat_Intensity = zeros(N_slices_total,7);
    feat_Hist = zeros(N_slices_total,6);
    feat_SVD = zeros(N_slices_total,780);
    feat_GLCM = zeros(N_slices_total,672);
    feat_RunLength = zeros(N_slices_total,44);
    feat_Fractal = zeros(N_slices_total,27);
    feat_FormFactor = zeros(N_slices_total,32);
    feat_Fourier = zeros(N_slices_total,300);
    feat_DCT = zeros(N_slices_total,2901);
    feat_Hankel = zeros(N_slices_total,75);
    feat_DistTrafo = zeros(N_slices_total,56);
    feat_TopHat = zeros(N_slices_total,60);
    feat_Skeleton = zeros(N_slices_total,17);
    feat_Unitary = zeros(N_slices_total,73);
    if strcmp(houghtype,'both')
        feat_Hough = zeros(N_slices_total,393);
    end
    if strcmp(houghtype,'linear')
        feat_Hough = zeros(N_slices_total,381);
    end
    if strcmp(houghtype,'circular')
        feat_Hough = zeros(N_slices_total,12);
    end
    feat_Zernike = zeros(N_slices_total,92);
    feat_Hu = zeros(N_slices_total,8);
    feat_Affine = zeros(N_slices_total,6);
    feat_LBP = zeros(N_slices_total,1024);
    feat_MSER = zeros(N_slices_total,15);
    feat_SalientRegion = zeros(N_slices_total,8);
    feat_Quadtree = zeros(N_slices_total,5);
    feat_EBR_IBR = zeros(N_slices_total,1211);
    feat_Connectivity = zeros(N_slices_total,47);
    feat_Harris = zeros(N_slices_total,10);
    feat_LineProfile = zeros(N_slices_total,122);
    feat_Law = zeros(N_slices_total,58);
    feat_LoG = zeros(N_slices_total,261);
    feat_Gilles = zeros(N_slices_total,6);
    feat_SURF = zeros(N_slices_total,11);
    feat_LOSIB = zeros(N_slices_total,34);
    feat_RCovD = zeros(N_slices_total,15);
    feat_Sector = zeros(N_slices_total,5);
elseif typeflag.global
    % if only all global features should be extracted
    feat_Intensity = zeros(N_slices_total,7);
    feat_Hist = zeros(N_slices_total,6);
    feat_SVD = zeros(N_slices_total,780);
    feat_GLCM = zeros(N_slices_total,672);
    feat_RunLength = zeros(N_slices_total,44);
    feat_Fractal = zeros(N_slices_total,27);
    feat_FormFactor = zeros(N_slices_total,32);
    feat_Fourier = zeros(N_slices_total,300);
    feat_DCT = zeros(N_slices_total,2901);
    feat_Hankel = zeros(N_slices_total,75);
    feat_DistTrafo = zeros(N_slices_total,56);
    feat_TopHat = zeros(N_slices_total,60);
    feat_Skeleton = zeros(N_slices_total,17);
    feat_Unitary = zeros(N_slices_total,73);
    if strcmp(houghtype,'both')
        feat_Hough = zeros(N_slices_total,393);
    end
    if strcmp(houghtype,'linear')
        feat_Hough = zeros(N_slices_total,381);
    end
    if strcmp(houghtype,'circular')
        feat_Hough = zeros(N_slices_total,12);
    end
    feat_Zernike = zeros(N_slices_total,92);
    feat_Hu = zeros(N_slices_total,8);
    feat_Affine = zeros(N_slices_total,6);
elseif typeflag.local
    % if only all local features should be extracted
    feat_LBP = zeros(N_slices_total,1024);
    feat_MSER = zeros(N_slices_total,15);
    feat_SalientRegion = zeros(N_slices_total,8);
    feat_Quadtree = zeros(N_slices_total,5);
    feat_EBR_IBR = zeros(N_slices_total,1211);
    feat_Connectivity = zeros(N_slices_total,47);
    feat_Harris = zeros(N_slices_total,10);
    feat_LineProfile = zeros(N_slices_total,122);
    feat_Law = zeros(N_slices_total,58);
    feat_LoG = zeros(N_slices_total,261);
    feat_Gilles = zeros(N_slices_total,6);        
    feat_Sector = zeros(N_slices_total,5);    
end


%% Preprocess images and extract features

iCounter = 1;

for iI = 1:length(images)
    I_3D = images{iI};
    %     image3D = Struct_Final.TestSample;
    
    for iSlice = 1:size(I_3D,3)
        
        % preprocessing, if necessary
        if doSegmentation || doGrayscaling
            for j = 1:size(I_3D,3)
                % Segmentation
                if doSegmentation
                    SegmentedImage = Segmentation(I_3D(:,:,j), 10 , 1000, 4,0);
                end
                % Convert image to 255 gray scale image
                if doGrayscaling
                    I = floor(mat2gray(SegmentedImage)*255);
                end
            end
        end
        
        I = I_3D(:,:,iSlice);
        
        %% feature extraction
        
        % #################################################################
        % GLOBAL FEATURES
        % #################################################################
        
        
        % -----------------------------------------------------------------
        % Intensity features
        % -----------------------------------------------------------------
        
        % Intensity-based features
        if (typeflag.global || typeflag.texture || typeflag.corr || typeflag.entropy)
            feat_Intensity(iCounter,:) = IntensityF(I,typeflag);
        end
        
        % Histogram-based features
        if (typeflag.global || typeflag.texture || typeflag.entropy)
            feat_Hist(iCounter,:) = HistogramF(I,typeflag);
        end
        
        % Singular Value Decompostion
        if (typeflag.global || typeflag.texture)
            feat_SVD(iCounter,:) = SVDF(I);
        end
        

        % -----------------------------------------------------------------
        % Geometrical features
        % -----------------------------------------------------------------
        
        % Gray Level Co-occurance Matrix
        if (typeflag.global || typeflag.texture || typeflag.corr || typeflag.entropy)
            feat_GLCM(iCounter,:) = GLCMF(I, GLCMParameters, typeflag);
        end
        
        % Run-Length
        if (typeflag.global || typeflag.texture)
            feat_RunLength(iCounter,:) = RunLengthF(I);
        end
        
        % Fractal Dimensions
        if (typeflag.global || typeflag.form)
            feat_Fractal(iCounter,:) = FractalDimensionF(I,plotflag);
        end

        % Form Factor
        if (typeflag.global || typeflag.form || typeflag.corr)
            feat_FormFactor(iCounter,:) = FormFactorF(I,typeflag);
        end
        
        % -----------------------------------------------------------------
        % Transformation features
        % -----------------------------------------------------------------
               
        % feat_Fouier
        if (typeflag.global || typeflag.transform || typeflag.corr || typeflag.moments)
            feat_Fourier(iCounter,:) = FourierTrafoF(I,typeflag);
        end
        
        % DCT
        if (typeflag.global || typeflag.transform || typeflag.corr)
            feat_DCT(iCounter,:) = DCTF(I,typeflag);
        end
        
        % Hankel
        if (typeflag.global || typeflag.transform || typeflag.corr)
            feat_Hankel(iCounter,:) = DHankelF(I,typeflag);
        end
        
        % Distance Transform
        if (typeflag.global || typeflag.transform || typeflag.corr)
            feat_DistTrafo(iCounter,:) = DistanceTrafoF(I,typeflag);
        end
        
        
        % Top-hat Transform
        if (typeflag.global || typeflag.transform || typeflag.form || typeflag.moments)
            feat_TopHat(iCounter,:) = TopHatTrafoF(I,typeflag);
        end
        
        % Skeletonization
        if (typeflag.global || typeflag.transform || typeflag.corr)
            feat_Skeleton(iCounter,:) = SkeletonizationF(I,typeflag);
        end
        
        % various other transforms
        if (typeflag.global || typeflag.transform)
            feat_Unitary(iCounter,:) = UnitaryTrafoF(I);
        end
        
        % Hough Transform
        if (typeflag.global || typeflag.transform || typeflag.form ||...
                (typeflag.moments && ~strcmp(houghtype,'circular')))
            feat_Hough(iCounter,:) = HoughTrafoF(I,houghtype,arc_min,plotflag,typeflag);
        end    
        
        % -----------------------------------------------------------------
        % Moment features
        % -----------------------------------------------------------------
        
        if (typeflag.global || typeflag.moments)
            % Zernike
            feat_Zernike(iCounter,:) = ZernikeF(I);
            
            % Hu
            feat_Hu(iCounter,:) = HuF(I);
            
            % Affine moment
            feat_Affine(iCounter,:) = AffineMomentsF(I);
        end
        
        
        % #################################################################
        % LOCAL FEATURES
        % #################################################################
        
        % -----------------------------------------------------------------
        % Region features
        % -----------------------------------------------------------------
        
        % Local Binary Pattern
        if (typeflag.local || typeflag.texture)
            feat_LBP(iCounter,:) = LocalBinaryPatternF(I, RadialLBP);
        end
        
        % MSER
        if (typeflag.local || typeflag.texture)
            feat_MSER(iCounter,:) = MSERF(I,0);
        end
        
        % Salient Region
        if (typeflag.local || typeflag.texture || typeflag.moments)
            feat_SalientRegion(iCounter,:) = SalientRegionF(I,typeflag);
        end
        
        % Quadtree Decomposition
        if (typeflag.local || typeflag.texture)
            feat_Quadtree(iCounter,:) = QuadtreeDecompositionF(I,threshold_quadtree);
        end
        
        % EBR and IBR
        if (typeflag.local || typeflag.form || typeflag.texture || typeflag.corr)
            feat_EBR_IBR(iCounter,:) = EBRandIBRF(I,typeflag);
        end
        
        
        % Connectivity
        if typeflag.local || typeflag.texture
            feat_Connectivity(iCounter,:) = ConnectivityF(I);
        end

        % Sector Decomposition
        if typeflag.local
            feat_Sector(iCounter,:) = SectorF(I);
        end
        
        % -----------------------------------------------------------------
        % Line features
        % -----------------------------------------------------------------
        
        % Harris Detector
        if (typeflag.local || typeflag.texture)
            feat_Harris(iCounter,:) = HarrisF(I,plotflag);
        end
        
        % Line profile
        if (typeflag.local || typeflag.texture || typeflag.moments || typeflag.corr)
            feat_LineProfile(iCounter,:) = LineProfileF(I,plotflag, typeflag);
        end
        
        % -----------------------------------------------------------------
        % Point features
        % -----------------------------------------------------------------
        
        % LAW
        if (typeflag.local || typeflag.texture || typeflag.moments)
            feat_Law(iCounter,:) = LawF(I);
        end
        
        % Laplacian of Gaussian
        if (typeflag.local || typeflag.texture || typeflag.moments)
            feat_LoG(iCounter,:) = LoGF(I,N_blobs,typeflag);
        end
        
        % Gilles Points
        if (typeflag.local || typeflag.entropy)
            feat_Gilles(iCounter,:) = GillesF(I);
        end
        
        % #################################################################
        % FEATURE DESCRIPTORS
        % #################################################################
        
        % SURF
        if typeflag.texture;
            feat_SURF(iCounter,:) = SURF(I);
        end
        
        % LOSIB
        if typeflag.texture
            feat_LOSIB(iCounter,:) = LOSIBF(I);
        end
        
        % Region of Covariance Descriptor
        if (typeflag.corr || typeflag.moments)
            feat_RCovD(iCounter,:) = RCovDsF(I,typeflag);
        end

        iCounter = iCounter + 1;
        
        
    end
end

%% Return feature vector for further image processing

if ~exist('feat_Intensity','var')
    feat_Intensity = [];
end    
if ~exist('feat_Hist','var')
    feat_Hist = [];
end    
if ~exist('feat_SVD','var')
    feat_SVD = [];
end    
if ~exist('feat_GLCM','var')
    feat_GLCM = [];
end    
if ~exist('feat_RunLength','var')
    feat_RunLength = [];
end    
if ~exist('feat_Fractal','var')
    feat_Fractal= [];
end    
if ~exist('feat_FormFactor','var')
    feat_FormFactor = [];
end    
if ~exist('feat_Fourier','var')
    feat_Fourier = [];
end    
if ~exist('feat_DCT','var')
    feat_DCT = [];
end    
if ~exist('feat_Hankel','var')
    feat_Hankel = [];
end    
if ~exist('feat_DistTrafo','var')
    feat_DistTrafo = [];
end    
if ~exist('feat_TopHat','var')
    feat_TopHat = [];
end    
if ~exist('feat_Skeleton','var')
    feat_Skeleton = [];
end    
if ~exist('feat_Unitary','var')
    feat_Unitary = [];
end    
if ~exist('feat_Hough','var')
    feat_Hough = [];
end    
if ~exist('feat_Zernike','var')
    feat_Zernike = [];
end    
if ~exist('feat_Hu','var')
    feat_Hu = [];
end    
if ~exist('feat_Affine','var')
    feat_Affine = [];
end    
if ~exist('feat_LBP','var')
    feat_LBP = [];
end    
if ~exist('feat_MSER','var')
    feat_MSER = [];
end    
if ~exist('feat_SalientRegion','var')
    feat_SalientRegion = [];
end    
if ~exist('feat_Quadtree','var')
    feat_Quadtree = [];
end    
if ~exist('feat_EBR_IBR','var')
    feat_EBR_IBR = [];
end    
if ~exist('feat_Connectivity','var')
    feat_Connectivity = [];
end    
if ~exist('feat_Harris','var')
    feat_Harris = [];
end    
if ~exist('feat_LineProfile','var')
    feat_LineProfile = [];
end    
if ~exist('feat_Law','var')
    feat_Law = [];
end    
if ~exist('feat_LoG','var')
    feat_LoG = [];
end    
if ~exist('feat_Gilles','var')
    feat_Gilles = [];
end    
if ~exist('feat_SURF','var')
    feat_SURF = [];
end    
if ~exist('feat_LOSIB','var')
    feat_LOSIB = [];
end    
if ~exist('feat_RCovD','var')
    feat_RCovD = [];
end    
if ~exist('feat_Sector','var')
    feat_Sector = [];
end    


feat_vector = [feat_Affine feat_RCovD feat_Connectivity feat_DCT...
    feat_DistTrafo feat_EBR_IBR feat_FormFactor feat_Fourier...
    feat_Gilles feat_Hankel feat_Harris feat_Hu feat_Law feat_LOSIB... 
    feat_LineProfile feat_LoG feat_MSER feat_Quadtree feat_SURF feat_SVD...
    feat_SalientRegion feat_TopHat feat_Unitary feat_Zernike feat_GLCM... 
    feat_Fractal feat_Hist feat_Intensity feat_LBP...
    feat_RunLength feat_Skeleton feat_Sector feat_Hough];


