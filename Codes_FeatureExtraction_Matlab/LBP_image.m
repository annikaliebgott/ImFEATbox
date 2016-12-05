function LBP = LBP_image(I, varargin) 
% Input:     - I: A 2D image
%            - filtR: a 2D matrix representing a round/radial filter. It
%              can be generated using the generateRadialFilterLBP function.
%            - isRotInv : a logical flag to enable rotation invariance in
%              the LBP image
%            - isChanWiseRot: a logical flag to enable channel wise 
%              rotation. 
%
% Output:     - LBP: LBP image
% 
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
%
%                           - The material by Nikolay S. in: 
%                           http://www.mathworks.com/matlabcentral/fileexchange/36484-local-binary-patterns


%% Deafult parameters
isRotInv = false;
isChanWiseRot = false;
filtR = generateRadialFilterLBP(8, 1);

%% Get user inputs overriding default values
funcParamsNames = {'filtR', 'isRotInv', 'isChanWiseRot'};
assignUserInputs(funcParamsNames, varargin{:});

% In case of file name input read graphical file
if ischar(I) && exist(I, 'file') == 2 
    I = imread(I);
end

nClrChans = size(I, 3);

inImgType = class(I);
calcClass = 'single';

isCalcClassInput = strcmpi(inImgType, calcClass);
if ~isCalcClassInput
    I = cast(I, calcClass);
end
imgSize = size(I);

nNeigh = size(filtR, 3);

if nNeigh <= 8
    outClass = 'uint8';
elseif nNeigh > 8 && nNeigh <= 16
    outClass = 'uint16';
elseif nNeigh > 16 && nNeigh <= 32
    outClass = 'uint32';
elseif nNeigh > 32 && nNeigh <= 64
    outClass = 'uint64';
else
    outClass = calcClass;
end

if isRotInv
    nRotLBP = nNeigh;
    nPixelsSingleChan = imgSize(1)*imgSize(2);
    iSingleChan = reshape(1:nPixelsSingleChan, imgSize(1), imgSize(2));
else
    nRotLBP = 1;
end

nEps = -3;
weigthVec = reshape(2.^((1:nNeigh)-1), 1, 1, nNeigh);
weigthMat = repmat(weigthVec, imgSize([1, 2]));
binaryWord = zeros(imgSize(1), imgSize(2), nNeigh, calcClass);
LBP = zeros(imgSize, outClass);
possibleLBP = zeros(imgSize(1), imgSize(2), nRotLBP);

for iChan = 1:nClrChans  
    % Initiate neighbours relation filter and LBP's matrix
    for iFiltElem = 1:nNeigh
        % Rotate filter to compare center to next neigbour
        filtNeight = filtR(:, :, iFiltElem);
        
        % calculate relevant LBP elements via filtering
        binaryWord(:, :, iFiltElem) = cast( ...
            roundn(filter2(filtNeight, I(:, :, iChan), 'same'), nEps) >= 0,...
            calcClass);
        % Without rounding sometimes inaqulity happens in some pixels
        % compared to pixelwiseLBP
    end 

    for iRot = 1:nRotLBP
        % find all relevant LBP candidates
        possibleLBP(:, :, iRot) = sum(binaryWord.*weigthMat, 3);
        
        % shift binaryWord elements
        if iRot < nRotLBP
            binaryWord = circshift(binaryWord, [0, 0, 1]); 
        end
    end
    
    if isRotInv
        if iChan == 1 || isChanWiseRot
            % Find minimal LBP, and the rotation applied to first color channel
            [minColroInvLBP, iMin] = min(possibleLBP, [], 3);
            
            % calculte 3D matrix index
            iCircShiftMinLBP = iSingleChan + (iMin-1)*nPixelsSingleChan;
        else
            % the above rotation of the first channel holds to rest of the channels
            minColroInvLBP = possibleLBP(iCircShiftMinLBP);
        end 
    else
        minColroInvLBP = possibleLBP;
    end 
    
    if strcmpi(outClass, calcClass)
        LBP(:, :, iChan) = minColroInvLBP;
    else
        LBP(:, :, iChan) = cast(minColroInvLBP, outClass);
    end
end

end