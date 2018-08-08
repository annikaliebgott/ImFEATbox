function img = fLoadImages(sPath, sFileName)
% function to load different image formats into a Matlab cell array
% compatible to ImFEATbox
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic
% and Interventional Radiology, University Hospital of Tuebingen, Germany
% and the Institute of Signal Processing and System Theory University of
% Stuttgart, Germany. Last modified: August 2018
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
%          thomas.kuestner@iss.uni-stuttgart.de
% ************************************************************************

if ~exist(sFileName)
    sFileName = [];
end

% catch error if there is no folder separator at the end of sPath
if ispc
    sep = '\';
else
    sep = '/';
end

if sPath(end)~= sep
    sPath = [sPath,sep];
end

% check if single file (e.g. one .jpg image) should be loaded or all files
% in one folder (e.g. all files of a 3D dicom image)
if isempty(sFileName)
    singleFile = false;
else
    singleFile = true;
end

% determine file type
if singleFile
    [~, ~, sType] = fileparts([sPath,sFileName]);
else
    files = dir(sPath);
    isDir = true;
    n = 1;
    while isDir
        isDir = files(n).isdir;
        n = n + 1;
    end
    [~, ~, sType] = fileparts([sPath,files(n-1).name]);
end

% load images
if singleFile
    if isempty(sType)
        disp('Warning: no file extension detected. Please check the files you want to load.')
        return
    else
        switch lower(sType)
            case {'.jpg', '.jpeg', '.png', '.bmp', '.tif', '.tiff'}
                img{1} = imread([sPath, sFileName]);
            case '.mat'
                data = load([sPath, sFileName]);
                datafield = fieldnames(data);
                img{1} = getfield(data,datafield{1});  
            case '.mhd'
                disp('to be added soon, sorry!')
                return
            case '.gipl'
                [giplImage,~] = fReadGIPL(sPath);
                img = {giplImage};    
            otherwise
                disp('Warning: file extension not supported. Please choose another file with one of the following formats:')
                disp('.jpg, .JPG, .jpeg, .JPEG, .bmp, .BMP, .tif, .TIF, .tiff, .TIFF, .mat, .dcm, .DCM, .nii, .NII, .mhd, .MHD, .gipl, .GIPL')
                return
        end
    end
    
else
    N_images = length(files)-sum([files.isdir]);
    N_dirs = sum([files.isdir]);
     
    if isempty(sType)
        disp('Warning: no file extension detected. Please check the files you want to load.')
        return
    else
        switch lower(sType)
            case {'.jpg', '.jpeg', '.png', '.bmp', '.tif', '.tiff'}
                img = cell(N_images,1);
                for m = N_dirs + 1 : length(files)
                    sFileName = files(m).name;
                    img{m-N_dirs} = imread([sPath, sFileName]);
                end
            case '.mat'
                img = cell(N_images,1);
                for m = N_dirs + 1 : length(files)
                    sFileName = files(m).name;
                    data = load([sPath, sFileName]);
                    datafield = fieldnames(data);
                    img{m-N_dirs} = getfield(data,datafield{1});
                end
            case {'.dcm','.ima'}
                [dicomImage,~,~] = fReadDICOM(sPath);
                img = {dicomImage};
            case '.nii'
                [niftyImage,~] = fReadNifty(sPath);
                img = {niftyImage};
                
            case '.mhd'
                disp('to be added soon, sorry!')
                return
            case '.gipl'
                [giplImage,~] = fReadGIPL(sPath);
                img = {giplImage};
            otherwise
                disp('Warning: file extension not supported. Please choose another file with one of the following formats:')
                disp('.jpg, .JPG, .jpeg, .JPEG, .bmp, .BMP, .tif, .TIF, .tiff, .TIFF, .mat, .dcm, .DCM, .nii, .NII, .mhd, .MHD, .gipl, .GIPL')
                return
        end
    end
    
    
end


end