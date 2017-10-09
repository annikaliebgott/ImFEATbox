# -*- coding: utf-8 -*-

import numpy as np
import skimage
from ImFEATbox.__helperCommands import grayscale2bw
from ImFEATbox.LocalFeatures.Region import findLocalMaximum

def GillesF(I, radius=10, th=0.95, returnShape=False):
    """
 Input:     - I: A 2D image
            - radius: radius of the region mask. Default: radius = 10
            - th: only keep points above a threshold defined by th.
              Default: th = 0.95

 Output:    - Out: A (1x6) vector containing 6 metrics calculated from
                   detected Gilles points
    """
# ************************************************************************
# Implemented for MRI feature extraction by the Department of Diagnostic
# and Interventional Radiology, University Hospital of Tübingen, Germany
# and the Institute of Signal Processing and System Theory University of
# Stuttgart, Germany. Last modified: November 2016
#
# Contact: annika.liebgott@iss.uni-stuttgart.de
# ************************************************************************

    if returnShape:
        return (6,1)

    if np.iscomplexobj(I):
        I = np.real(I)

    ## extract Gilles points

    # convert image to binary
    im = grayscale2bw(I)

    # define a region mask
    mask = skimage.morphology.disk(radius)

    # compute local entropy
    loc_entropy = skimage.filter.rank.entropy(im,mask)

    # find the local maxima
    #local_max = findLocalMaximum(loc_entropy,radius)[2]
    local_max = scipy.ndimage.maximum_filter(val, size=10, mode='constant', footprint=mask)
    # keep only points above a threshold
    row,col = np.where(local_max > th*np.max(local_max))

    # normalize coordinates for better comparison of different sized images
    points_gilles = np.array([row, col])/float(I.size)

    ## feature extraction
    points_gilles_num = np.shape(points_gilles)[0]
    if points_gilles_num > 1:
        points_gilles_mean = np.mean(points_gilles, axis=0)
        points_gilles_std = np.std(points_gilles, axis=0, ddof=1)
        points_gilles_std2 = np.std(points_gilles, ddof=1)
    else:
        if points_gilles_num == 0:
            points_gilles_mean = np.array([0, 0])
        else:
            points_gilles_mean = points_gilles
        points_gilles_std = np.array([0, 0])
        points_gilles_std2 = 0

    ## return features
    Out = np.hstack([points_gilles_num, points_gilles_mean, points_gilles_std, points_gilles_std2])

    return Out
