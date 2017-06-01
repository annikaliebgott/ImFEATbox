function [SegmentedImage, SegmentedBackground, BinaryImage] = Segmentation(I, Margin, cSegParam, iDim, lCrop, plotflag)
% Function to segment the MR images based on Chan-Vese model
%
% Input     - I: A 2D/3D image. Color images will be converted to grayscale                
%           - Margin: Margin values from the edges of the image for 
%             rectangular initial masking. Binary image representing 
%             initialization. Initial state position being close to object.
%           - cParam: Chan-Vese parameters:
%                     - cParam{1} = Lambda 1: inside
%                     - cParam{2} = Lambda 2: outside
%                     - cParam{3} = mu: step size
%                     - cParam{4} = Smoothness: Smoothnes Factor. Higher value leads to smoother
%                                   boundaries, smooth out details
%                     - cParam{5} = Iteration: Number of iterarions to run for Chan-Vese model                      
%           - lCrop: crop image to a reduced size (with zero-padding around)
%           - plotflag: flag to enable/disable visualization with to be plotted slice
%
% Output    - SegmentedImage: The main image with black backgrounds
%           - SegmentedBackground
%           - BinaryImage: The black and white map of segmentation.
%            (1: foreground, 0: background)
%
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic
% and Interventional Radiology, University Hospital of Tuebingen, Germany
% and the Institute of Signal Processing and System Theory University of
% Stuttgart, Germany. Last modified: June 2017
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
% Modified: thomas.kuestner@iss.uni-stuttgart.de
% ************************************************************************
%
% Implementation based on:  T. F. Chan and L. A. Vese, "Active contours 
%                           without edges", Image processing, 
%                           IEEE transactions on, vol. 10, no. 2, 
%                           pp. 266-277, 2001."

if(nargin < 6), plotflag = false; end;
if(nargin < 5), lCrop = true; end;
if(nargin < 4), iDim = ndims(I); end;

% Check for color image and convert to the grayscale
if(numel(size(I))==3 && iDim ~= 3)       %if image is 3D
    if(size(I,3)==3)
        I = rgb2gray(I);
    end
end

% Initialization
MaskImage = zeros(size(I));
MaskImage(Margin:size(I,1)-Margin,Margin:size(I,2)-Margin)=1;
MaskImage = double(MaskImage);


% Check whether the MaskImage and original image sizes are the same and
% raise an error
if(any(size(I)~=size(MaskImage)))
    error('Image and MaskImage Must Be in The Same Size');
end

if(iDim == 2) % 2D segmentation
    % Perform segmentation by means of Chan-Vese model MATLAB function.
    BinaryImage = activecontour(abs(I),MaskImage, cSegParam{5}, 'Chan-Vese',cSegParam{4});
elseif(iDim == 3) % 3D segmentation
    BinaryImage = ChanVese('exe', abs(I), MaskImage, [], cSegParam);
else
    error('Segmentation(): Undefined segmentation algorithm');
end

if(lCrop)
    [ SegmentedImage, SegmentedBackground ] = fSegmentCrop ( I, BinaryImage );
else % keep image size
    SegmentedImage(BinaryImage) = I(BinaryImage);
    SegmentedBackground(~BinaryImage) = I(~BinaryImage);
end

% plotflag = false;
if plotflag > 0
    % display results
    figure; subplot(2,2,1)
    imagesc(I(:,:,plotflag)); axis image; colormap gray;
    title('The original image');
    
    subplot(2,2,2)
    imagesc(MaskImage(:,:,plotflag)); axis image; colormap gray;
    title('The initialization mask');
    
    subplot(2,2,3)
    imagesc(BinaryImage(:,:,plotflag)); axis image; colormap gray;
    title({'The final', 'segmenatation output'});
    
    subplot(2,2,4)
    imagesc(I(:,:,plotflag)); axis image; colormap gray;
    hold on;
    contour(BinaryImage(:,:,plotflag),'b','linewidth',3);
    hold off;
    title({'The image with the', 'segmentation (blue)'});
end
