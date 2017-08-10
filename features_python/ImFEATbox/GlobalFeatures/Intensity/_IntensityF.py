# -*- coding: utf-8 -*-

import numpy as np
from ImFEATbox.__helperCommands import rgb2grayscale

def IntensityF(I, typeflag=None, returnShape=False):
    """
     Input:     - I: A 2D image
                - typeflag: dict of logicals to permit extracting features
                  based on desired characteristics:
                       + typeflag['global']: all features
                       + typeflag['texture']: all features
                       + typeflag['corr']: only features based on correlation
                       + typeflag['entropy']: only features based on entropy
                  default: all features are being extracted
                  For more information see README.txt

     Output:    - Out: A (1x7) vector containing 7 metrics based on intensity
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
    #
    # Implementation based on the intensity-based features of the paper by
    # McGee et al.: "Image metric-based correction (Autocorrection) of motion
    # effects: Analysis of image metrics" (Journal of Magnetic Resonance
    # Imaging, vol. 11, no. 2, pp. 174–181, 2000.)

    if typeflag == None:
        typeflag = dict()
        typeflag['global'] = True
        typeflag['texture'] = True
        typeflag['corr'] = True
        typeflag['entropy'] = True

    if typeflag['global'] or typeflag['texture']:
        typeflag['corr'] = True
        typeflag['entropy'] = True

    # Check for color image and convert to grayscale
    if len(np.shape(I)) == 3:
        if np.shape(I)[2] == 3:
            I = rgb2grayscale(I)

    Height, Width = np.shape(I)
    n = Height * Width

    ## Extract features

    if typeflag['global'] or typeflag['texture'] or typeflag['corr']:
        # AutoCorrelation 1
        p2 = np.zeros(Height,Width-1)
        p1 = np.sum(np.power(I,2))
        k = range(0, Width)
        for i in range(0, Height):
            p2[i,k] = I[i,k]*I[i,k+1]
        ACORR = p1 - np.sum(p2[:])

        # AutoCorrelation 2
        p3 = np.zeros(Height,Width-1)
        l = range(0, Width-1)
        for i in range(0, Height):
            p3[i,l] = I[i,l]*I[i,l+1]
        ACORR2 = np.sum(p2) - np.sum(p3)

    if typeflag['global'] or typeflag['texture'] or typeflag['entropy']:
        # Entropy
        H1 = 1.0*I / (np.sum(np.power(I, 2)))
        # marginal entropies
        E = -np.sum(H1[H1 != 0] * np.log(H1[H1 != 0]))

    if typeflag['global'] or typeflag['texture']:

        #Standard Deviation
        STD = np.std(I)

        # Cube of Normalized Intensities
        NI3 = np.sum(np.power((1.0*I)/np.sum(I),3))

        # 4th power of Normalized Intensities
        NI4 = np.sum(np.power((1.0*I)/np.sum(I),4))

        # Squared Intensities
        I2 = np.sum(np.power((1.0*I)/n,2))

    ## return feature vector
    if not (typeflag['texture'] or typeflag['global']):
        if not typeflag['corr']:
            Out = E
        elif not typeflag['entropy']:
            Out = [ACORR, ACORR2]
        else:
            Out = [ACORR, ACORR2, E]
    else:
        Out = [STD, ACORR, ACORR2, E, NI3, NI4, I2]

    return Out
