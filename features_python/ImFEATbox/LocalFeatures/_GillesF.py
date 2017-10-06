import numpy as np
import skimage

def GillesF(I, radius=10, th=0.95):
    """
# Input:     - I: A 2D image
#            - radius: radius of the region mask. Default: radius = 10
#            - th: only keep points above a threshold defined by th.
#              Default: th = 0.95
#
# Output:    - Out: A (1x6) vector containing 6 metrics calculated from
#                   detected Gilles points
    """
# ************************************************************************
# Implemented for MRI feature extraction by the Department of Diagnostic
# and Interventional Radiology, University Hospital of Tübingen, Germany
# and the Institute of Signal Processing and System Theory University of
# Stuttgart, Germany. Last modified: November 2016
#
# Contact: annika.liebgott@iss.uni-stuttgart.de
# ************************************************************************



    if np.iscomplexobj(I):
        I = np.real(I)

    ## extract Gilles points

    # convert image to binary
    im = I < skimage.filter.threshold_otsu(I)

    # define a region mask
    mask = skimage.morphology.disk(r)

    # compute local entropy
    loc_entropy = skimage.filter.rank.entropy(im,mask)

    # find the local maxima
    [~,~,local_max] = findLocalMaximum(loc_entropy,radius)

    # keep only points above a threshold
    [row,col] = find(local_max > th*max(local_max(:)))

    # normalize coordinates for better comparison of different sized images
    points_gilles = [row col]/numel(I)

    ## feature extraction
    points_gilles_num = size(points_gilles,1)
    if points_gilles_num > 1
        points_gilles_mean = mean(points_gilles)
        points_gilles_std = std(points_gilles)
        points_gilles_std2 = std2(points_gilles)
    else
        if points_gilles_num == 0
            points_gilles_mean = [0 0]
        else
            points_gilles_mean = points_gilles
        end
        points_gilles_std = [0 0]
        points_gilles_std2 = 0
    end

    ## return features
    Out = [points_gilles_num points_gilles_mean points_gilles_std points_gilles_std2]

    end
