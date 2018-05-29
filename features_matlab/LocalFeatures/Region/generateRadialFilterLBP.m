function [radInterpFilt] = generateRadialFilterLBP(p, r)
% Input:       - p: an integer number specifying the number of neigbours, 
%                number of enabled filter elements
%              - r: a positive number specifying the raduis of the filter
%
%
% Output:     - radInterpFilt: filter to generate LBP image with function
%               LBP_image.m
% 
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
% Implementation based on:  - T. Ojala, M. Pietikainen and T. Maenpaa, 
%                           "Multiresolution gray-scale and rotation 
%                           invariant texture classification with local 
%                           binary patterns", Pattern Analysis and Machine
%                           Intelligence, IEEE Transactions on, vol. 24, 
%                           no. 7, pp. 971–987, 2002.
%
% Implemented by:       Nikolay S. 
%                       First version: Nikolay S. 2014-01-09.
%                       Last update:   Nikolay S. 2014-01-16.
%                       see also: http://www.mathworks.com/matlabcentral/fileexchange/36484-local-binary-patterns


%% Example
% [radInterpFilt]=generateRadialFilterLBP(40, 20);
% figure;
% imshow( sum(radInterpFilt, 3)~=0 );
% title('Wighted filter');
%

%% Default parameters
if nargin<2
    r=1;
    if nargin<1
        p=8;
    end
end

%% verify parameters legit values
% radius and number of neighbours must not be below 1, number of neighbors
% must be an integer
r=max(1, r);   
p=round(p);     
p=max(1, p);   


%% find elements angles, aranged counter clockwise starting from "X axis"
% See http://www.ee.oulu.fi/mvg/files/pdf/pdf_6.pdf for illustration

theta=linspace(0, 2*pi, p+1)+pi/2;   

% remove obsolete last element (0=2*pi)
theta=theta(1:end-1);           


%% Find relevant coordinates

% convert to cartesian
[rowsFilt, colsFilt] = pol2cart(theta, repmat(r, size(theta))); 
nEps=-3;
rowsFilt=roundn(rowsFilt, nEps);
colsFilt=roundn(colsFilt, nEps);

% Matrix indices should be integers
rowsFloor=floor(rowsFilt);
rowsCeil=ceil(rowsFilt);
colsFloor=floor(colsFilt);
colsCeil=ceil(colsFilt);

rowsDistFloor=1-abs(rowsFloor-rowsFilt);
rowsDistCeil=1-abs(rowsCeil-rowsFilt);
colsDistFloor=1-abs(colsFloor-colsFilt);
colsDistCeil=1-abs(colsCeil-colsFilt);

% Find minimal filter dimentions, based on indexes
filtDims=[ceil(max(rowsFilt))-floor(min(rowsFilt)),...
    ceil(max(colsFilt))-floor(min(colsFilt))];
% verify filter dimentions are odd
filtDims=filtDims+mod(filtDims+1, 2);

filtCenter=(filtDims+1)/2;


%% Convert cartesian coordinates to matrix element coordinates via simple shift
rowsFloor=rowsFloor+filtCenter(1);
rowsCeil=rowsCeil+filtCenter(1);
colsFloor=colsFloor+filtCenter(2);
colsCeil=colsCeil+filtCenter(2);


%% Generate the filter, each 2D slice for filter element  
radInterpFilt=zeros([filtDims,  p], 'single'); 
for iP=1:p
    radInterpFilt(rowsFloor(iP), colsFloor(iP), iP)=...
        radInterpFilt(rowsFloor(iP), colsFloor(iP), iP)+rowsDistFloor(iP)+colsDistFloor(iP);
    
    radInterpFilt(rowsFloor(iP), colsCeil(iP), iP)=...
        radInterpFilt(rowsFloor(iP), colsCeil(iP), iP)+rowsDistFloor(iP)+colsDistCeil(iP);
    
    radInterpFilt(rowsCeil(iP), colsFloor(iP), iP)=...
        radInterpFilt(rowsCeil(iP), colsFloor(iP), iP)+rowsDistCeil(iP)+colsDistFloor(iP);
   
    radInterpFilt(rowsCeil(iP), colsCeil(iP), iP)=...
        radInterpFilt(rowsCeil(iP), colsCeil(iP), iP)+rowsDistCeil(iP)+colsDistCeil(iP);
    
    radInterpFilt(:, :, iP)=radInterpFilt(:, :, iP)/sum(sum(radInterpFilt(:, :, iP)));
end

% Substract 1 at central element to get difference between central element and relevant
% neighbours: (5) T=p{s(g1-g0), s(g2-g0),...,s(gn-g0)}
radInterpFilt(filtCenter(1), filtCenter(2), :)=...
    radInterpFilt(filtCenter(1), filtCenter(2), :)-1; 

end