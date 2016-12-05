function Out = LocalBinaryPatternF(I, RadialInput)
% Input:     - I: A 2D image
%            - RadialInput: A (nx2) array containing the radius and number 
%              of pixels to generate the LBP circular filter. n is the  
%              number of filters, the first column is the number of pixels  
%              and the second column is the radius.
%              Example : RadialInput = [8 1; 16 2]
%
%
% Output:    - Out: A (1xn*256) vector containing n*256 LBP histograms from 
%                   filtering the image with different filters 
% 
% Note: This implementation calls the functions LBP_image.m and
%       generateRadialFilterLBP.m, which are also included in ImFEATbox.
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
% Implementation based on:  - T. Ojala, M. Pietikainen and T. Maenpaa, 
%                           "Multiresolution gray-scale and rotation 
%                           invariant texture classification with local 
%                           binary patterns", Pattern Analysis and Machine
%                           Intelligence, IEEE Transactions on, vol. 24, 
%                           no. 7, pp. 971–987, 2002.


m = 1;
for n =1:size(RadialInput,1)
    LBPArray = LBP_image(I, 'filtR', generateRadialFilterLBP(RadialInput(n,1), RadialInput(n,2)));
    LBPHist (m:m+256-1) = histc(LBPArray(:),0:255)';
    m = 256 + m; 
end

Out = LBPHist;

end