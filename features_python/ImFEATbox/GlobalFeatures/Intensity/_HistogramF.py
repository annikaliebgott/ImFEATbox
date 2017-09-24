# -*- coding: utf-8 -*-

import numpy as np
from ImFEATbox.__helperCommands import rgb2grayscale

def HistogramF(I, typeflag=None, returnShape=False):
    """
     Input:     - I: A 2D image
                - typeflag: Struct of logicals to permit extracting features
                  based on desired characteristics:
                       + typeflag['global']: all features
                       + typeflag['texture']: all features
                       + typeflag['entropy']: only features based on entropy
                  default: all features are being extracted
                  For more information see README.txt

     Output:    - Out: A (1x6) vector containing 6 metrics calculated from the
                       image histogram
   """
    # ************************************************************************
    # Implemented for MRI feature extraction by the Department of Diagnostic
    # and Interventional Radiology, University Hospital of Tuebingen, Germany
    # and the Institute of Signal Processing and System Theory University of
    # Stuttgart, Germany. Last modified: November 2016
    #
    # This implementation is part of ImFEATbox, a toolbox for image feature
    # extraction and analysis. Available online at:
    # https://github.com/annikaliebgott/ImFEATbox
    #
    # Contact: annika.liebgott@iss.uni-stuttgart.de
    # ************************************************************************

    if typeflag == None:
        typeflag = dict()
        typeflag['global'] = True
        typeflag['texture'] = True
        typeflag['entropy'] = True

    if typeflag['global'] or typeflag['texture']:
        typeflag['entropy'] = True

    if returnShape:
        if typeflag['texture'] or typeflag['global']:
            return (6,1)
        else:
            return (1,1)


    # Check for color image and convert to grayscale
    if len(np.shape(I)) == 3:
        if np.shape(I)[2] == 3:
            I = rgb2grayscale(I)

    # 256 gray scale Image
    graylevels = range(0, 256)

    # Probability of occurence of gray values
    Prob = np.histogram(I, bins=len(graylevels), range=(np.min(0),np.max(I)))[0] / float(I.size)

    ## extract features

    if typeflag['texture'] or typeflag['global']:
        # Histogram Mean _ Identifier code
        Mean = np.dot(graylevels, Prob)

        # Histogram Standard Deviation
        Std = np.sqrt(np.dot(np.power((graylevels-Mean),2),Prob))

        # Histogram Skewness
        Skewness = 1/np.power(Std,3)*np.dot(np.power((graylevels-Mean),3),Prob)

        # Histogram Kurtosis
        Kurtosis = 1/np.power(Std,4)*(np.dot(np.power((graylevels-Mean),4),Prob)-3)

        # Histogram Energy
        Energy = np.sum(np.power(Prob,2))

    # Histogram Entropy
    # marginal entropies
    Entropy = -np.sum(Prob[Prob != 0] * np.log2(Prob[Prob != 0]))


    ## Return feature vector

    if typeflag['texture'] or typeflag['global']:
        Out = np.array([Mean, Std, Skewness, Kurtosis, Energy, Entropy])
    else:
        Out =  np.array([Entropy])

    return Out
