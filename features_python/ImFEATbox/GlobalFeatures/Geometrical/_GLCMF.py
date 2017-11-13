# -*- coding: utf-8 -*-

import numpy as np
from ImFEATbox.__helperCommands import rgb2grayscale, isColorImage
from skimage.feature import greycomatrix, greycoprops
from skimage import exposure
from scipy.stats import entropy

def GLCMF(I, typeflag=None, DisplacementVector=np.array([1]), NumLevels=8, GrayLimits=None, returnShape=False):
    """
        Input:     - I: A 2D image
                   - DisplacementVector: A (nx2) vector composed of offset
                       and orientation (default = [0 1])
                   - NumLevels: An integer number specifying the number of
                       gray level intensities (default = 8)
                   - GrayLimits: A (1x2) vector containing the min/max
                       gray levels which are needed to sort gray levels of I
                       into number of gray levels specified by NumLevels
                   - typeflag: Struct of logicals to permit extracting features
                        based on desired characteristics:
                       + typeflag['global']: all features
                       + typeflag['form']: all features
                       + typeflag['texture']: all features
                       + typeflag['corr']: only features based on correlation
                       + typeflag['entropy']: only features based on entropy
                  default: all features are being extracted
                  For more information see README.txt

     Output:    - Out: A (1xn*21) vector containing n*21 metrics calculated
                  from the gray level co-occurence matrix
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
# Implementation based on:  R. Haralick; K. Shanmugam; I. Dinstein (1973):
#                           "Textural Features for Image Classification".
#                           IEEE Transactions on Systems, Man, and
#                           Cybernetics. SMC-3 (6): 610–621.



    if GrayLimits == None:
        GrayLimits = np.array([np.min(I), np.max(I)])


    # Angles could be implemented in the future as input parameter
    Angles=None
    if Angles == None:
        Angles = np.array([0])

    I = np.array(I)

    # Check for color image and convert to grayscale

    if isColorImage(I):
        I = rgb2grayscale(I)

    # graycomatrix.m can't process complex input values
    #Image = double(real(Image));
    if np.any(np.iscomplex(I)):
        I = I.real

    # Default Parameters (Offset value = 1 pixel)

#    if 'typeflag' not in globals():
    # if ~exist('InputParameters','var') or ~isfield(InputParameters,'DisplacementVector')
    #     InputParameters.DisplacementVector = [0 1];
    # end

    # if ~exist('InputParameters','var') or ~isfield(InputParameters,'NumLevels')
    #     InputParameters.NumLevels = 8;
    # end
    #
    # if ~exist('InputParameters','var') or ~isfield(InputParameters,'GrayLimits')
    #     InputParameters.GrayLimits = [min(Image(:)) max(Image(:))];
    # end

    if typeflag == None:
        typeflag = dict()
        typeflag['global'] = True
        typeflag['texture'] = True
        typeflag['form'] = True
        typeflag['corr'] = True
        typeflag['entropy'] = True


    if typeflag['global'] or typeflag['texture']:
        typeflag['corr'] = True
        typeflag['entropy'] = True


    numDist = len(DisplacementVector)

    if returnShape:
        if typeflag['texture'] or typeflag['global']:
            return (numDist * 21, 1)
        elif typeflag['corr']:
            if typeflag['entropy']:
                return (numDist * 7, 1)
            else:
                return (numDist * 4, 1)
        else:
            return (numDist * 5, 1)






    uI = exposure.equalize_hist(I, nbins=8)
    uI = uI-uI.min()
    uI = 7*uI/float(np.max(uI))
    uI = uI.astype('uint8')
    DisplacementVector = np.array([1])
    #DisplacementVector = np.array([1])
    #Angles = [0, np.pi/2]
    Angles = np.array([0])
    #NumLevels = 8



    # Image must be uint8
    GLCM_Matrices = greycomatrix(image = uI,
        distances = DisplacementVector,
        angles = Angles,
        levels = NumLevels,
        # TODO limits wichtig??
        #'GrayLimits', InputParameters.GrayLimits,
        symmetric = False)
        # levels x levels x number of distances x number of angles.


    GLCM_Matrices = GLCM_Matrices[:,:,:,0] # angles are not implemented yet


    nMax = np.shape(GLCM_Matrices)[2]

    if typeflag['texture'] or typeflag['global']:
        # extract all features of GLCM
        Out = np.zeros((nMax, 21))
    elif typeflag['corr']:
        if typeflag['entropy']:
            # ectract correlation based and entropy based features
            Out = np.zeros((nMax, 7))
        else:
            # extract only correlation based features
            Out = np.zeros((nMax, 4))
    else:
        # extract only entropy based features
        Out = np.zeros((nMax, 5))


    #for n=1:1:size(GLCM_Matrices,3)
    for n in range(nMax):

        G = GLCM_Matrices[:,:,n]
        SumOfGLCM = np.sum(G)
        GNormalized = G / float(SumOfGLCM)
        MeanOfGLCM = np.mean(GNormalized)
        s1, s2 = np.shape(G)[0:2]

        # Initialization
        p_x = np.zeros(s1)
        p_y = np.zeros(s2)
        u_x = 0
        u_y = 0
        HXY1 = 0
        HXY2 = 0
        sigma_x = 0
        sigma_y = 0
        p_xplusy = np.zeros(s1*2 - 1)
        p_xminusy = np.zeros(s1)
        CO = 0
        CORR = 0
        IDM = 0
        SE = 0
        H = 0
        DE = 0
        IMC1 = 0
        IMC2 = 0
        ACORR = 0
        DIS = 0
        CLS = 0
        CLP = 0
        INV = 0
        INVN = 0
        IDMN = 0


        m = np.array(range(1, s1+1)) #replace one for-loop by using a vector of length s1
        o = np.array(range(1, s2+1)) #replace one for-loop by using a vector of length s2

        # Angular Second Moment or Energy (ASM)
        ASM = np.sum(np.power(GNormalized, 2))

        # Sum of squared variance (SSV)
        SSV = np.sum(np.power(m - MeanOfGLCM, 2) * GNormalized)

        # Entropy (H)
        if typeflag['entropy']:
            #H = - np.sum(GNormalized.flatten() * np.log(GNormalized.flatten()))
            H = entropy(GNormalized.flatten())


        #for i = 1:s1
        for i in range(0, s1):
            # TODO hier weiter! i korrekt anpassen!
            u_x = u_x + np.sum((i+1)*GNormalized[i,:])
            u_y = u_y + np.sum(o * GNormalized[i,:])

            p_x[i] = np.sum(GNormalized[i,:])
            p_y[i] = np.sum(GNormalized[:,i])

            #k = (i+1)+(1:s2)
            #l = np.abs(i-(1:s2))
            p_xplusy[i:i+s2] = p_xplusy[i:i+s2] + GNormalized[i,:s2].T

            # since there are multiple equal values in index vector l (i.e. some
            # indices are in there twice) and matlab can't automatically assign
            # the sum of multiple values to the same index, the calculation has
            # to be divided into two parts (on from first index down until l==1,
            # then up from l==0 to l(end)
            print(i)
            p_xminusy[:i][::-1] = p_xminusy[:i][::-1] + GNormalized[i,:i]
            p_xminusy[:s2-i] = p_xminusy[:s2-i] + GNormalized[i,i:s2]

            # Contrast (CO)
            CO = CO + np.sum(np.power((i+1 - o), 2) * GNormalized[i,:])

            # Inverse difference moment or homogenity (IDM)
            IDM = IDM + np.sum(GNormalized[i,o-1] / ( 1 + np.power(i+1 - o.astype('float'),2)))

            # Dissimilarity (DIS)
            DIS = DIS + np.sum(np.abs(i+1 - o) * GNormalized[i,:])

            # Autocorrelation (ACORR)
            if typeflag['corr']:
                ACORR = ACORR + np.sum((i+1) * o * GNormalized[i,o])

            # Inverse difference (INV)
            INV = INV + np.sum(GNormalized[i,o-1] / ( 1 + np.abs(i+1 - o.astype('float'))))

            # Inverse difference normalized (INVN)
            INVN = INVN + np.sum(GNormalized[i,o-1] / ( 1 + (np.abs(i+1-o) / float(s1))))

            # Inverse difference moment normalized (IDMN)
            IDMN = IDMN + np.sum(GNormalized[i,o-1] / ( 1 + np.power(i+1-o / float(s1),2)))


        #for i = 1:s1
        for i in range(s1):

            # Cluster shade (CLS)
            CLS = CLS + (np.power(i + o - u_x - u_y,3)) * GNormalized[i,o-1].T

            # Cluster prominence (CLP)
            CLP = CLP + (np.power(i + o - u_x - u_y,4)) * GNormalized[i,o-1].T

            sigma_x = sigma_x  + np.sum((np.power((i+1) - u_x, 2)) * GNormalized[i,o-1])
            sigma_y = sigma_y  + np.sum((np.power((o) - u_y, 2)) * GNormalized[i,o-1])


        # Correlation (CORR)
        if typeflag['corr']:
            CORR = (ACORR - (u_x*u_y)) / float(np.sqrt(sigma_x)*np.sqrt(sigma_y))


        # Summed average (SA)
        #SA = (2:2*s1)*p_xplusy
        SA = range(2,2*s1+1) * p_xplusy

        # Summed entropy (SE)
        if typeflag['entropy']:
            SE = -np.sum(p_xplusy * np.log(p_xplusy))

        # Summed Variance (SV)
        SV = (np.power(np.array(range(2,2*s1+1)) - SE, 2)) * p_xplusy

        # Difference varience (DV)
        DV = np.power(range(0,s1), 2) * p_xminusy

        # Difference entropy (DE)
        if typeflag['entropy']:
            DE = - np.sum(p_xminusy * np.log(p_xminusy))

            HXY = H

            for i in range(s1):
                HXY1 = HXY1 - GNormalized[i,:] * np.log(p_x[i] * p_y)
                HXY2 = HXY2 - np.sum(p_x[i] * p_y * np.log(p_x[i] * p_y))

            Hx = - np.sum(p_x * np.log(p_x))
            Hy = - np.sum(p_y * np.log(p_y))

            # Information measure of correlation 1 (IMC1)
            # Information measure of correlation 2 (IMC2)
            IMC1 = (HXY - HXY1) / float(np.max([Hx, Hy]))
            IMC2 = np.power(1 - np.exp(-2*(HXY2 - HXY)), 0.5)

        # Maximum probability (MAXP)
        MAXP = np.max(GNormalized)

        if typeflag['texture'] or typeflag['global']:
            Out[n,:] = np.array([ACORR , CO, CORR, CLP, CLS ,  DIS, ASM, H, IDM, MAXP,
                SSV, SA, SV, SE, DV, DE, IMC1, IMC2, INV,INVN, IDMN])
        elif typeflag['corr']:
            if typeflag['entropy']:
                Out[n,:] = np.array([ACORR, CORR, H, SE, DE, IMC1, IMC2])
            else:
                Out[n,:] = np.array([ACORR, CORR, IMC1, IMC2])
        else:
            Out[n,:] = np.array([H, SE, DE, IMC1, IMC2])

    print(Out)

    return np.hstack(Out)
