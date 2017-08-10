%% segmentation
% 2D/3D segmentation (2D = false, 3D = true)
lSegDim = false;

% perform segmentation and subsequently process:
% 0: no segmentation -> process whole image
% 1: segmentation -> process foreground image
% 2: segmentation -> process background image
doSegmentation = 1;

% Chan-Vese parameters
% lambda1, lambda2, mu, smoothness, iterations
cSegParam = {1, 1, 0.2, 4, 100};
% margin from edges
iMargin = 10;

%% scaling
% scale image into gray scale value range 
% doGrayscaling=0: no scaling
% scalar doGraysacling ~= 0: scaling into range [0 doGraysacling]
% vector doGraysacling: scaling into range [doGraysacling(1) doGraysacling(2)]
% doGrayscaling<0: 2D scaling (slice-wise)
% doGrayscaling>0: 3D scaling (volume)
doGrayscaling = 255;

%% feature parameters
% GLCM
values = [1;2;4;8;16;32;64;128];
offset0 = [zeros(size(values,1),1) values];
offset1 = [-values values];
offset2 = [-values zeros(size(values,1),1)];
offset3 = [-values -values];
offset = [offset0 ; offset1 ; offset2 ; offset3];
GLCMParameters.DisplacementVector = offset;
GLCMParameters.NumLevels = 255;

% LocalBinaryPatternF.m: set filter parameters for LBP feature extraction
RadialLBP = [8 1; 16 2;24 3;32 4];

% LoGF.m: set desired number of blobs for Laplacian of Gaussian feature
% extraction
N_blobs = 120;

% QuadtreeDecompositionF.m: set threshold for Quadtree decomposition
% feature extraction
threshold_quadtree = 0.27;

% HoughTrafoF.m: set type and minimum arc for Hough transform features
houghtype = 'both';
arc_min = pi/2;

% HarrisF.m: choose number of points considered to be the strongest points
N_s_Harris = 25;

% FractalDimensionF.m: choose largest box size
width = 512;

% GillesF.m: define mask radius and threshold factor
radius_Gilles = 10;
threshold_Gilles = 0.95;

% LoGF.m: define LoG-parameters begin, end and step size for sigma_array
sigma.begin = 2;
sigma.end = 15;
sigma.step = 1;

% SURF.m: choose number of points considered to be the strongest points
N_s_SURF = 25;

% SalientRegionF.m: set border of rectangles for Haar feature calculation
% and inner Region
border_Salient = 14;
R1 = 1; 

% TopHatTrafoF.m: create morphological structruning element (STREL)
SE = cell(zeros);
SE{1} = strel('arbitrary', [1 0 0;1 0 0;1 0 1]);
SE{2} = strel('ball', 15,5);
SE{3} = strel('diamond', 13);
SE{4} = strel('disk', 15);
SE{5} = strel('line', 10, 45);
SE{6} = strel('octagon', 3);
SE{7} = strel('pair', [2 2]);
SE{8} = strel('periodicline', 2, [1 -2]);
SE{9} = strel('rectangle', [3 5]);
% SE{10} = strel('square', 11);

% UnitaryTrafoF.m: choose which transformations you wish to carry out
transformation.WH = false;
transformation.H = false;
transformation.R = false;
transformation.C = true;

% GradientF.m: specify which gradient (first and/or second order gradients)
% you wish to apply
gradtype.first = true;
gradtype.second = true;

% GaborFilterF.m: set parameters for Gabor filters
scale = 5;
orientation = 8;

% LacunarityF.m: set parameters for box counting alogorithm of lacunarity
l_min = 10;
l_max = 30;
stepsize = 1;

% WaveletTrafoF.m: specify which type of wavelet you wish to employ
wavelettype.haar = true;
wavelettype.coif = true;
wavelettype.dmey = true;

% LOSIBF.m: specify desired radius and number of neighbors
radius = 5;
neighbors = 16;

%% debugging
% set whether visualizations included in some feature extraction algorithms
% should be plotted or not (note: increases computation time)
plotflag = false;