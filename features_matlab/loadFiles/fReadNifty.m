function [dImg, dDim] = fReadNifty(sFilename)
% ************************************************************************
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
%          thomas.kuestner@iss.uni-stuttgart.de
% ************************************************************************
% Implemented by and Copyright: 
% Copyright 2012-2015 Christian Wuerslin, Stanford University: 
% IMAGINE: IMAGe visualization, analysis and evaluation engINE
%   Code adapted from: "Read Medical Data 3D" by Dirk-Jan Kroon (c) 2010:
%   https://www.mathworks.com/matlabcentral/fileexchange/29344-read-medical-data-3d

fid = fopen(sFilename, 'rb');
if(fid < 0)
    error('Could not open the file ''%s''!\n', sFilename);
end

% -------------------------------------------------------------------------
% Read the relevant hearder data
fseek(fid, 42, 'bof');
iSize = fread(fid, 7, 'uint16');
iSize(iSize == 0) = 1;
fseek(fid, 14, 'cof');
iType = fread(fid, 1, 'uint16');
fseek(fid, 8, 'cof');
dDim  = fread(fid, 7, 'float');
fseek(fid, 4, 'cof');
dSlope = fread(fid, 1, 'float');
dIntercept = fread(fid, 1, 'float');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Find start of volume data
iNumEl = prod(iSize);
switch iType
    case 1,                     dBytesPerVoxel = 1./8;
    case {2, 256},              dBytesPerVoxel = 1;
    case {4, 512},              dBytesPerVoxel = 2;
    case {128},                 dBytesPerVoxel = 3;
    case {8, 16, 768, 2304},    dBytesPerVoxel = 4;
    case {32, 64, 1024, 1280},  dBytesPerVoxel = 8;
    case {1536, 1792},          dBytesPerVoxel = 16;
    case {2048},                dBytesPerVoxel = 32;
    otherwise,                  error('Unknown datatype');
end
fseek(fid, -iNumEl.*dBytesPerVoxel, 'eof');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Read volume data (convert all data to double)
switch iType
    case     1, dImg = fread(fid, iNumEl, 'bit1');
    case     2, dImg = fread(fid, iNumEl, 'uint8');
    case     4, dImg = fread(fid, iNumEl, 'int16');
    case     8, dImg = fread(fid, iNumEl, 'int32');
    case    16, dImg = fread(fid, iNumEl, 'float');
    case    32, dImg = fread(fid, iNumEl.*2, 'float'); %Complex
    case    64, dImg = fread(fid, iNumEl, 'double');
    case   128, dImg = fread(fid, iNumEl.*3, 'uint8'); % RGB
    case   256, dImg = fread(fid, iNumEl, 'int8');
    case   512, dImg = fread(fid, iNumEl, 'uint16');
    case   768, dImg = fread(fid, iNumEl, 'uint32');
    case  1024, dImg = fread(fid, iNumEl, 'int64');
    case  1280, dImg = fread(fid, iNumEl, 'uint64');
    case  1792, dImg = fread(fid, iNumEl.*2, 'double'); %Complex
    otherwise, error('Format not suported');
end
fclose(fid);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Get the image into the right format
dImg = dImg.*dSlope + dIntercept;
switch iType
    case {32, 1792} % Complex
        dImg = complex(dImg(1:2:end), dImg(2:2:end));
    case 128 %RGB
        warning('RGB Mode is not supported. Converting to grayscale!');
        dImg = reshape(dImg, [3, length(dImg)./3]);
        dImg = mean(dImg);
end
dImg = reshape(dImg(:), iSize');
dDim = dDim(1:ndims(dImg));
% -------------------------------------------------------------------------%