function img = fLoadImages(sPath, sFileName)
% Function to load different image formats into a Matlab cell array
% compatible to ImFEATbox. The following file formats are currently
% supported:
% Category 1: .jpg, jpeg, .png, .tif, .tiff, .bmp
% Category 2: .mat
% Category 3: .nii
% Category 4: .mhd
% Category 5: .gipl
% Category 6: .dcm
%
% 1.) To load single files, please use img = fLoadImages(sPath, sFileName) 
% with sPath = path to folder containing the file and sFileName = name of 
% the file you wish to load. 
%
% 2.) To load multiple files at once, use img = fLoadImages(sPath,[]).
% Note: right now, it is only possible to load multiple files if there are 
% other files in the folder (i.e. no unsupported files like .pdf) and all
% files to be loaded belong to the same image category.
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
            case '.nii'
                [niftyImage,~] = fReadNifty([sPath, sFileName]);
                img = {niftyImage};   
            case '.mhd'
                [data,~] = fReadMHD([sPath, sFileName]);
                datafield = fieldnames(data);
                img{1} = getfield(data, datafield{1});
            case '.gipl'
                [giplImage,~] = fReadGIPL([sPath, sFileName]);
                img = {giplImage};    
            otherwise
                disp('Warning: file extension not supported. Please choose another file with one of the following formats:')
                disp('.jpg, .JPG, .jpeg, .JPEG, .bmp, .BMP, .tif, .TIF, .tiff, .TIFF, .mat, .dcm, .DCM, .nii, .NII, .mhd, .MHD, .gipl, .GIPL')
                return
        end
    end
    
else
    N_dirs = sum([files.isdir]);
     
    if isempty(sType)
        disp('Warning: no file extension detected. Please check the files you want to load.')
        return
    else
        switch lower(sType)
            case {'.jpg', '.jpeg', '.png', '.bmp', '.tif', '.tiff'}
                N_images = length(files)-sum([files.isdir]);
                img = cell(N_images,1);
                for m = N_dirs + 1 : length(files)
                    sFileName = files(m).name;
                    img{m-N_dirs} = imread([sPath, sFileName]);
                end
            case '.mat'
                N_images = length(files)-sum([files.isdir]);
                img = cell(N_images,1);
                for m = N_dirs + 1 : length(files)
                    sFileName = files(m).name;
                    data = load([sPath, sFileName]);
                    datafield = fieldnames(data);
                    img{m-N_dirs} = getfield(data,datafield{1});
                end
            case '.nii'
                N_images = length(files)-sum([files.isdir]);
                img = cell(N_images,1);
                for m = N_dirs + 1 : length(files)
                    sFileName = files(m).name;
                    [niftyImage,~] = fReadNifty([sPath,sFileName]);
                    img{m-N_dirs} = niftyImage;      
                end    
            case '.mhd'
                N_images = (length(files)-sum([files.isdir]))/2;
                img = cell(N_images,1);
                k = 1;
                for m = N_dirs + 1 : length(files)               
                    sFileName = files(m).name;
                    [~, ~, sTypeFile] = fileparts([sPath,sFileName]);
                    if strcmpi(sTypeFile,'.mhd')
                        [data,~] = fReadMHD([sPath, sFileName]);
                        datafield = fieldnames(data);
                        img{k} = getfield(data, datafield{1});
                        k = k + 1;
                    end    
                end
            case '.gipl'
                N_images = length(files)-sum([files.isdir]);
                img = cell(N_images,1);
                for m = N_dirs + 1 : length(files)
                    sFileName = files(m).name;
                    [giplImage,~] = fReadGIPL([sPath, sFileName]);
                    img{m-N_dirs} = giplImage;
                end    
            case {'.dcm','.ima'}
                [dicomImage,~,~] = fReadDICOM(sPath);
                img = {dicomImage};    
            otherwise
                disp('Warning: file extension not supported. Please choose another file with one of the following formats:')
                disp('.jpg, .JPG, .jpeg, .JPEG, .bmp, .BMP, .tif, .TIF, .tiff, .TIFF, .mat, .dcm, .DCM, .nii, .NII, .mhd, .MHD, .gipl, .GIPL')
                return
        end
    end
    
    
end


end