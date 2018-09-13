function Out = LOSIBF(I, radius, neighbors)
% Input:     - I: A 2D image
%            - radius: radius to consider around a pixel.
%              default: radius = 5                
%            - neighbors: Number of neighbors to consider
%              default: neighbors = 16
%
% Output:    - Out: A (1x34) vector containing 34 metrics 
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
% Implemented by:   Oscar Garcia-Olalla (ogaro@unileon.es)
%                   VARP Group (http://pitia.unileon.es/VARP)
%                   University of Leon
%                   Contact email: ogaro@unileon.es
%                   Modified:27/11/14 Revision: 1.0
%
% Reference:    O. García-Olalla, E. Alegre, L. Fernández-Robles and 
%               V. González-Castro, "Local Oriented Statistics Information 
%               Booster (LOSIB) for Texture Classification," Pattern 
%               Recognition (ICPR), 2014 22nd International Conference on, 
%               Stockholm, 2014, pp. 1114-1119.



% Default parameters
if ~exist('radius','var')
    radius = 5;
end
if ~exist('neighbors','var')
    neighbors = 16; 
end

%% LOSIB

% Initializate the variable that holds the neighbours coordinates.
neighborPoints = zeros(neighbors,2);

% Angle between each neighbor.
angle = 2*pi/neighbors; 

% Determine the coordinates of the neighbors given their angle and radius.
for i = 1:neighbors
    neighborPoints(i,1) = -radius*sin((i-1)*angle);
    neighborPoints(i,2) = radius*cos((i-1)*angle);
end

% Convert color image to gray scale
if (size(I,3)==3)
    I = rgb2gray(I);
end
I = im2double(I);
[numRows,numCols] = size(I);

numRowsPatch = radius*2 + 1; % --radius--C--radius--
numColsPatch = radius*2 + 1; % --radius--C--radius--

% Determine the first valid pixel (does not exceed the limits of the image).
originRow = 1 + radius;
originCol = 1 + radius; 

% Number of rows/columns computed (all except the border)
dRows = numRows - numRowsPatch; 
dCols = numCols - numColsPatch;

% store valid center pixels of the image
imageCenters = I(originRow:originRow+dRows, originCol:originCol+dCols);

% Preallocate array to store neighbor point images in
NArray = zeros(numel(imageCenters),neighbors+1);

for i = 1:neighbors
    row = neighborPoints(i,1)+originRow;
    col = neighborPoints(i,2)+originCol;
    
    % Calculate floors, ceils and rounds for the x and y.
    floorRow = floor(row); 
    ceilRow = ceil(row); 
    roundRow = round(row);
    floorCol = floor(col); 
    ceilCol = ceil(col); 
    roundCol = round(col);
    
    % Check if interpolation is needed
    if (abs(row - roundRow) == 0) && (abs(col - roundCol) ==0)
        % Interpolation is not needed, use original datatypes
        % Select a matrix with upper left corner the pixel of the neighbour
        % and dimension dRows+1 X dCols+1
        N = I(roundRow:roundRow+dRows,roundCol:roundCol+dCols); 
        % convert the matrix to a column vector
        NArray(:,i) = N(:);
    else
        % Interpolation needed, use double type images
        
        % Get just the decimal part
        ty = row - floorRow; 
        tx = col - floorCol; 
        
        % Calculate the interpolation weights.
        w1 = (1 - tx) * (1 - ty);
        w2 =      tx  * (1 - ty);
        w3 = (1 - tx) *      ty ;
        w4 =      tx  *      ty ;
        
        % Compute interpolated pixel values. Sum up the matrix of the four 
        % points nearer the interpolated point multiplying them by w.
        N = w1*I(floorRow:floorRow+dRows,floorCol:floorCol+dCols)... 
            + w2*I(floorRow:floorRow+dRows,ceilCol:ceilCol+dCols)... 
            + w3*I(ceilRow:ceilRow+dRows,floorCol:floorCol+dCols)... 
            + w4*I(ceilRow:ceilRow+dRows,ceilCol:ceilCol+dCols);
        
        % convert the matrix to a column vector
        NArray(:,i) = N(:);
    end
end

% In the last column, we allocate the image with all the centers.
NArray(:,neighbors+1) = imageCenters(:);

losib_mean = zeros(1,neighbors);
losib_median = zeros(1,neighbors);

% Calculate the difference between the central pixel and the corresponding 
% neighbor for all the pixels in the image. Also, calculate mean and median
% of the difference vector as feature
for i=1:neighbors
    difference = abs(NArray(:,i)-NArray(:,neighbors+1));
    losib_mean(i) = mean(difference);
    losib_median(i) = median(difference);
end

% standard deviation
losib_mean_std = std(losib_mean);
losib_median_std = std(losib_median);

%% return feature vector
Out = [losib_mean losib_median losib_mean_std losib_median_std];

end
