def fExtractFeatures(image, lMask, cFeatures, sFeatureParas):
    """
         Main function to carry out the feature extraction

         This function calls all feature extraction algorithms provided by
         ImFEATbox.

         Input:    image: 2D/3D image to be processed
                   lMask: predefined masking
                   cFeatures: cell containing
                       a) the wanted features (with names provided by '-getFeatures')
                       b) complete feature groups: 'all' | 'global' | 'local' | 'corr' | 'gradient' | 'moments' | 'texture' | 'form' | 'entropy' | 'transform'
                       c) '-getFeatures': to retrieve implemented features
                   sFeatureParas: string with path to parameter file script; if not specified using default parameters
         Output:   feat_vec: an [NxM] array of M features extracted from N images.
                   feat_names: [2xM] cell array of M feature names and their column index

         The basic steps to successfully carry out the feature extraction:
         1.)   Import the images of which you want to extract the features.
               Note: this script expects images to be cell arrays
         2.)   Set typeflag to define which features you wish to extract.
               For more information about typeflag: see README.txt
         3.)   Set parameters needed for some feature extraction algorithms.
               Note: default values are tuned for use with automated MR image
               quality assessment, you might need to change them according to your
               application
         4.)   Choose wheather or not you wish to use additional tools
               (preprocessing, visualization).
    """
    # ************************************************************************
    # Implemented for MRI feature extraction by the Department of Diagnostic
    # and Interventional Radiology, University Hospital of Tuebingen, Germany
    # and the Institute of Signal Processing and System Theory University of
    # Stuttgart, Germany. Last modified: February 2017
    #
    # This implementation is part of ImFEATbox, a toolbox for image feature
    # extraction and analysis. Available online at:
    # https://github.com/annikaliebgott/ImFEATbox
    #
    # Contact: annika.liebgott@iss.uni-stuttgart.de, thomas.kuestner@iss.uni-stuttgart.de
    # ************************************************************************


    ## initialization
    currpath = fileparts(mfilename('fullpath'))
    addpath(genpath(currpath))

    if(nargin < 3 || ~exist('sFeatureParas','var') || ~ischar(sFeatureParas))
        sFeatureParas = [currpath,filesep,'parameters_ImFEATBox_def.m']
    end
    [pathPara, filenamePara, ~] = fileparts(sFeatureParas)
    cd(pathPara)
    eval([filenamePara,''])
    cd(currpath)
    if(ischar(cFeatures))
        cFeatures = {cFeatures}
    end

    # check if graphical elements can be shown
    if(usejava('jvm') && ~feature('ShowFigureWindows'))
        lDisplay = false
        plotflag = false
    else
        lDisplay = true
    end

    ## parse feature algos
    # we want:
    # hFeatures handleMethods, call functions of the featrues
    # cUsedFeatures: complete line of feature.. only using the name later


    [hFeatures, cUsedFeatures, typeflag] = fParseFeatureAlgos( cFeatures, sFeatureParas)
    if(isempty(hFeatures)) # shortcut
        feat_vector = cUsedFeatures
        feat_names = cUsedFeatures
        return
    end


    ## preprocessing
    # preallocate feature arrays for speed
    # for iI=1:length(hFeatures)
    #     eval(sprintf('feat_#s = zeros(size(image,3), #d)', cUsedFeatures{iI,2}, cUsedFeatures{iI,4}))
    # end

    image3D = image
    # segmentation
    if(isempty(lMask))
        if doSegmentation > 0
            if(lDisplay), multiWaitbar( 'Segmentation', 0 ) end
            if(lSegDim) # 3D
                [SegmentedImage, SegmentedBackground] = Segmentation(image3D, iMargin , cSegParam, 3, true, plotflag) # 3D segmentation
                if(lDisplay), multiWaitbar( 'Segmentation', 'Value', 1 ) end
            else # 2D
                BinaryImage = zeros(size(image3D))
                for iSlice = 1:size(image3D,3)
                    [~, ~, BinaryImage(:,:,iSlice)] = Segmentation(image3D(:,:,iSlice), iMargin , cSegParam, 2, false, plotflag) # 2D segmentation
                    if(lDisplay), multiWaitbar( 'Segmentation', 'Value', iSlice/size(image3D,3) ) end
                end
                [ SegmentedImage, SegmentedBackground ] = fSegmentCrop ( image3D, BinaryImage )
            end
            if(doSegmentation == 1)
                image3D = SegmentedImage
            elseif(doSegmentation == 2)
                image3D = SegmentedBackground
            end
            if(lDisplay), multiWaitbar( 'Segmentation', 'Close' ) end
        end
    else
        # precrop image according to mask
        image3D = fSegmentCrop ( image3D, lMask )
    end

    # scaling
    if(any(doGrayscaling ~= 0))
        # Scaling: convert image to N gray scale values
        if(any(doGrayscaling < 0)) # 2D (slice-wise) scaling
            if(isscalar(doGrayscaling))
                iRange = [0 -doGrayscaling]
            else
                iRange = -doGrayscaling
            end
            for iSli=1:size(image3D,3)
                image3D(:,:,iSli) = scaleImg(image3D(:,:,iSli),iRange)
            end
        else # 3D (volume) scaling
            if(isscalar(doGrayscaling))
                iRange = [0 doGrayscaling]
            else
                iRange = doGrayscaling
            end
            image3D = scaleImg(image3D,iRange)
        end
    end

    ## feature extraction
    # Display warning message in case I contains complex numbers
    if any(imag(image3D(:))~=0)
       disp('Warning: image contains complex numbers.')
       disp('Some feature extraction algorithms can not process complex values and therefore only use the real part.')
    end
    feat_vector = []
    feat_names = cell(2,length(hFeatures))
    # 2D feature extraction
    if(lDisplay), multiWaitbar( 'Slice', 0 ) end
    for iSlice = 1:size(image3D,3)
        I = image3D(:,:,iSlice)

        feat_sli = []
        for iI = 1:length(hFeatures)
            if(lDisplay && iSlice == 1), multiWaitbar(  cUsedFeatures{iI,1}, 0 ) end
            hFeat = hFeatures{iI}
            feat_curr = hFeat(I)
            if(iSlice == 1)
                feat_names{1,iI} = cUsedFeatures{iI,1}
                feat_names{2,iI} = [length(feat_sli)+1,length(feat_sli)+length(feat_curr)]
            end
            feat_sli = cat(2,feat_sli, feat_curr)
    #         eval(sprintf('feat_#s(iSlice,:) = hFeat(I)', cUsedFeatures{iI,2}))
            if(lDisplay), multiWaitbar( cUsedFeatures{iI,1}, 'Value', iSlice/size(image3D,3) ) end
        end
        feat_vector = cat(1,feat_vector,feat_sli)
        if(lDisplay), multiWaitbar( 'Slice', 'Value', iSlice/size(image3D,3) ) end
    end
    if(lDisplay), multiWaitbar( 'CloseAll' ) end

    return [feat_vector, feat_names]
