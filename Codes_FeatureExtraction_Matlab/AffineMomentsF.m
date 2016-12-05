function Out = AffineMomentsF(I)
% Input:    - I: A 2D image
%
%
% Output:   - Out: A (1x6) vector containing 6 moment features 
%
% ************************************************************************
% Modified for MRI feature extraction by the Department of Diagnostic 
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
% Implementation based on:  Tomas Suk, Jan Flusser, "Combined Blur and 
%                           Affine Moment Invariants and their use in 
%                           Pattern Recognition", Pattern Recognition, 
%                           vol. 36, 2003.
%
% Implemented by:   Asad Ali. Email: m.aliasad@yahoo.com



[x,y,pixelValues]=find(I(:,:,1));

pixelValues = double(pixelValues);
m00 = sum(pixelValues);

x = x - sum(x.*pixelValues)/m00;
y = y - sum(y.*pixelValues)/m00;

%% calculate moments

% second order central moments
m20 = CentralMoments(x,y,2,0,pixelValues);
m02 = CentralMoments(x,y,0,2,pixelValues);
m11 = CentralMoments(x,y,1,1,pixelValues);

% third order central moments
m30 = CentralMoments(x,y,3,0,pixelValues);
m03 = CentralMoments(x,y,0,3,pixelValues);
m21 = CentralMoments(x,y,2,1,pixelValues);
m12 = CentralMoments(x,y,1,2,pixelValues);

% fouth order central moments
m40 = CentralMoments(x,y,4,0,pixelValues);
m04 = CentralMoments(x,y,0,4,pixelValues);
m31 = CentralMoments(x,y,3,1,pixelValues);
m13 = CentralMoments(x,y,1,3,pixelValues);
m22 = CentralMoments(x,y,2,2,pixelValues);

% fifth order central moments
m50 = CentralMoments(x,y,5,0,pixelValues);
m05 = CentralMoments(x,y,0,5,pixelValues);
m41 = CentralMoments(x,y,4,1,pixelValues);
m14 = CentralMoments(x,y,1,4,pixelValues);
m32 = CentralMoments(x,y,3,2,pixelValues);
m23 = CentralMoments(x,y,2,3,pixelValues);

% seventh order central moments
m70 = CentralMoments(x,y,7,0,pixelValues);
m07 = CentralMoments(x,y,0,7,pixelValues);
m16 = CentralMoments(x,y,1,6,pixelValues);
m61 = CentralMoments(x,y,6,1,pixelValues);
m52 = CentralMoments(x,y,5,2,pixelValues);
m25 = CentralMoments(x,y,2,5,pixelValues);
m43 = CentralMoments(x,y,4,3,pixelValues);
m34 = CentralMoments(x,y,3,4,pixelValues);


% for blur invariance we recompute certain values
m50 = m50 - (10*m30*m20/m00);
m41 = m41 - (2*(3*m21*m20 + 2*m30*m11)/m00);
m32 = m32 - ((3*m12*m20 + m30*m02 + 6*m21*m11)/m00);
m23 = m23 - ((3*m21*m02 + m03*m20 + 6*m12*m11)/m00);
m14 = m14 - (2*(3*m12*m02 + 2*m03*m11)/m00);
m05 = m05 - (10*m03*m02/m00);

% for blur invariance seventh order moments recomputed
m70 = m70 - 7 * (3*m50*m20 + 5*m30*m40)/m00 + (210*m30*m20^2 / m00^2);
m61 = m61 - (6*m50*m11 + 15*m41*m20 + 15*m40*m21 + 20*m31*m30)/m00 + ...
    30*(3*m21*m20^2 + 4*m30*m20*m11)/m00^2;
m52 = m52 - (m50*m02 +10*m30*m22 + 10*m32*m20 + 20*m31*m21 +10*m41*m11 + 5*m40*m12)/m00 + ...
    10* (3*m12*m20^2 + 2*m30*m20*m02 + 4*m30*m11^2 + 12*m21*m20*m11)/m00^2;

m43 = m43 - (m40*m03 + 18*m21*m22 + 12*m31*m12 + 4*m30*m13 + 3*m41*m02 + 12*m32*m11 + ...
    6*m23*m20)/m00 + 6*(m03*m20^2 + 4*m30*m11*m02 + 12*m21*m11^2 + 12*m12*m20*m11 + 6*m21*m02*m20);

m34 = m34 - (m04*m30 + 18*m12*m22 + 12*m13*m21 + 4*m03*m31 + 3*m14*m20 + 12*m23*m11 ...
    + 6*m32*m02)/m00 + 6 *(m30*m02^2 + 4*m03*m11*m20 + 12*m12*m11^2 + 12*m21*m02*m11 + ...
    6*m12*m20*m02)/m00^2;

m25 = m25 - (m05*m20 + 10*m03*m22 + 10*m23*m02 + 20*m13*m12 + 10*m14*m11 + 5*m04*m21)/m00 + ...
    10*(3*m21*m02^2 + 2*m03*m02*m20 +4*m03*m11^2 + 12*m12*m02*m11)/m00^2;

m16 = m16 - (6*m05*m11 + 15*m14*m02 + 15*m04*m12 + 20*m13*m03)/m00 + 30*(3*m12*m02^2 + ...
    4*m03*m02*m11)/m00^2;

m07 = m07 - 7*(3*m05*m02 + 5*m03*m04)/m00 + (210*m03*m02^2 / m00^2); 

% first invariant computed from the determinant of the polynomial
I1 = (m30^2*m03^2 - 6*m30*m21*m12*m03 + 4*m30*m12^3 + ...
     4*m21^3*m03 - 3*m21^2*m12^2) / m00^10;

I2 = (m50^2*m05^2 - 10*m50*m41*m14*m05 + 4*m50*m32*m23*m05 + ...
    16*m50*m32*m14^2 - 12*m50*m23^2*m14 + 16*m41^2*m23*m05 + ...
    9*m41^2*m14^2 - 12*m41*m32^2*m05 - 76*m41*m32*m23*m14 + ...
    48*m41*m23^3 + 48*m32^3*m14 - 32*m32^2*m23^2)/m00^14;

I3 = (m30^2*m12*m05 - m30^2*m03*m14 - m30*m21^2*m05 - 2*m30*m21*m12*m14 + ...
    4*m30*m21*m03*m23 + 2*m30*m12^2*m23 - 4*m30*m12*m03*m32 + ...
    m30*m03^2*m41 + 3*m21^3*m14 - 6*m21^2*m12*m23 - 2*m21^2*m03*m32 + ...
    6*m21*m12^2*m32 + 2*m21*m12*m03*m41 - m21*m03^2*m50 - 3*m12^3*m41 + ...
    m12^2*m03*m50) / m00^11;

I4 = (2*m30*m12*m41*m05 - 8*m30*m12*m32*m14 + 6*m30*m12*m23^2 - ...
    m30*m03*m50*m05 + 3*m30*m03*m41*m14 - 2*m30*m03*m32*m23 - ...
    2*m21^2*m41*m05 + 8*m21^2*m32*m14 - 6*m21^2*m23^2 + ...
    m21*m12*m50*m05 - 3*m21*m12*m41*m14 + 2*m21*m12*m32*m23 + ...
    2*m21*m03*m50*m14 - 8*m21*m03*m41*m23 + 6*m21*m03*m32^2 - ...
    2*m12^2*m50*m14 + 8*m12^2*m41*m23 - 6*m12^2*m32^2)/m00^12;

I5 = (m30*m41*m23*m05 - m30*m41*m14^2 - m30*m32^2*m05 + 2*m30*m32*m23*m14 - ...
    m30*m23^3 - m21*m50*m23*m05 + m21*m50*m14^2 + m21*m41*m32*m05 - ...
    m21*m41*m23*m14 - m21*m32^2*m14 + m21*m32*m23^2 + m12*m50*m32*m05 - ...
    m12*m50*m23*m14 - m12*m41^2*m05 + m12*m41*m32*m14 + m12*m41*m23^2 - ...
    m12*m32^2*m23 - m03*m50*m32*m14 + m03*m50*m23^2 + ...
    m03*m41^2*m14 - 2*m03*m41*m32*m23 + m03*m32^3)/m00^13;

I6 = (m70^2*m07^2 - 14*m70*m61*m16*m07 + 18*m70*m52*m25*m07 + 24*m70*m52*m16^2 - ...
    10*m70*m43*m34*m07 - 60*m70*m43*m25*m16 + 40*m70*m34^2*m16 + 24*m61^2*m25*m07 + ...
    25*m61^2*m16^2 - 60*m61*m52*m34*m07 - 234*m61*m52*m25*m16 + 40*m61*m43^2*m07 + ...
    50*m61*m43*m34*m16 + 360*m61*m43*m25^2 - 240*m61*m34^2*m25 + 360*m52^2*m34*m16 + ...
    81*m52^2*m25^2 - 240*m52*m43^2*m16 - 990*m52*m43*m34*m25 + 600*m52*m34^3 + ...
    600*m43^3*m25 - 375*m43^2*m34^2)/m00^18;


%% return feature vector
Out = [I1 I2 I3 I4 I5 I6];

end

% Calculate Central Moments
function [cenMoment] = CentralMoments(x,y,p,q,pixelValues)
    cenMoment = sum(x.^p .* y.^q .* pixelValues);
end    
