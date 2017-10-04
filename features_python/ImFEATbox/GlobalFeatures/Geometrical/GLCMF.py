import numpy as np




def GLCMF(Image, DisplacementVector=np.array([0,1]), NumLevels=8, GrayLimits=np.array([np.min(Image), np.max(Image)]), typeflag):
""" Input:     - I: A 2D image
               - DisplacementVector: A (nx2) vector composed of offset
                   and orientation (default = [0 1])
               - NumLevels: An integer number specifying the number of
                   gray level intensities (default = 8)
               - GrayLimits: A (1x2) vector containing the min/max
                   gray levels which are needed to sort gray levels of I
                   into number of gray levels specified by NumLevels
               - typeflag: Struct of logicals to permit extracting features
                    based on desired characteristics:
                   + typeflag.global: all features
                   + typeflag.form: all features
                   + typeflag.texture: all features
                   + typeflag.corr: only features based on correlation
                   + typeflag.entropy: only features based on entropy
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


    Image = np.array(Image)

    # Check for color image and convert to grayscale

    if len(Image.shape) == 4 and Image.shape[3] == 3      #image is 3D and color
        if(size(Image,3)==3) # TODO FIX correct detection
            Image = rgb2gray(Image)





    # graycomatrix.m can't process complex input values
    #Image = double(real(Image));
    if np.any(np.iscomplex(Image)):
        Image = Image.real

    # Default Parameters (Offset value = 1 pixel)

#    if 'typeflag' not in globals():
    # if ~exist('InputParameters','var') || ~isfield(InputParameters,'DisplacementVector')
    #     InputParameters.DisplacementVector = [0 1];
    # end

    # if ~exist('InputParameters','var') || ~isfield(InputParameters,'NumLevels')
    #     InputParameters.NumLevels = 8;
    # end
    #
    # if ~exist('InputParameters','var') || ~isfield(InputParameters,'GrayLimits')
    #     InputParameters.GrayLimits = [min(Image(:)) max(Image(:))];
    # end

    if 'typeflag' not in globals():
    #if ~exist('typeflag', 'var')
        typeflag.global = True
        typeflag.texture = True
        typeflag.form = True
        typeflag.corr = True
        typeflag.entropy = True


    if typeflag.global || typeflag.texture:
        typeflag.corr = True
        typeflag.entropy = True


    # Image must be uint8
    GLCM_Matrices = graycomatrix(image = Image,
        distances = InputParameters.DisplacementVector
        levels = InputParameters.NumLevels,
        # TODO limits wichtig??
        #'GrayLimits', InputParameters.GrayLimits,
        symmetric = False)

    if typeflag.texture || typeflag.global
        # extract all features of GLCM
        Out = np.zeros(np.shape(GLCM_Matrices)[2], 21)
    elif typeflag.corr
        if typeflag.entropy
            # ectract correlation based and entropy based features
            Out = np.zeros(np.shape(GLCM_Matrices)[2], 7)
        else
            # extract only correlation based features
            Out = np.zeros(np.shape(GLCM_Matrices)[2], 4)
        end
    else
        # extract only entropy based features
        Out = np.zeros(np.shape(GLCM_Matrices)[2], 5)
    end


    #for n=1:1:size(GLCM_Matrices,3)
    for n in range(np.shape(GLCM_Matrices)[2]):

        G = GLCM_Matrices[:,:,n]
        SumOfGLCM = np.sum(G)
        GNormalized = G / float(SumOfGLCM)
        MeanOfGLCM = np.mean(GNormalized)
        s1, s2 = np.shape(G)[0:2]

        # Initialization
        p_x   = np.zeros(s1,1)
        p_y   = np.zeros(s2,1)
        u_x   = 0
        u_y   = 0
        HXY1  = 0
        HXY2  = 0
        sigma_x = 0
        sigma_y = 0
        p_xplusy = np.zeros((s1*2 - 1),1)
        p_xminusy = np.zeros(s1,1)
        CO      = 0
        CORR    = 0
        IDM     = 0
        SE      = 0
        H       = 0
        DE      = 0
        IMC1    = 0
        IMC2    = 0
        ACORR   = 0
        DIS     = 0
        CLS     = 0
        CLP     = 0
        INV     = 0
        INVN    = 0
        IDMN    = 0


        m = range(1, s1+1) #replace one for-loop by using a vector of length s1
        o = range(1, s2+1) #replace one for-loop by using a vector of length s2

        # Angular Second Moment or Energy (ASM)
        ASM = np.sum(np.power(GNormalized[:], 2))

        # Sum of squared variance (SSV)
        SSV = np.sum(np.power(m - MeanOfGLCM, 2) * GNormalized)

        # Entropy (H)
        if typeflag.entropy:
            H = - np.sum(np.dot(GNormalized[:], np.log(GNormalized[:])))
        end

        #for i = 1:s1
        for i in range(0, s1):
            # TODO hier weiter! i korrekt anpassen!
            u_x = u_x + np.sum((i+1)*GNormalized[i,:])
            u_y = u_y + np.sum(o * GNormalized[i,:])

            p_x[i] = np.sum(GNormalized[i,:])
            p_y[i] = np.sum(GNormalized[:,i])

            k = (i+1)+(1:s2)
            l = abs(i-(1:s2))
            p_xplusy[k-1] = p_xplusy[k-1] + GNormalized[i,(1:s2)].T

            # since there are multiple equal values in index vector l (i.e. some
            # indices are in there twice) and matlab can't automatically assign
            # the sum of multiple values to the same index, the calculation has
            # to be divided into two parts (on from first index down until l==1,
            # then up from l==0 to l(end)
            p_xminusy[l[1:i]+1] = p_xminusy[l[1:i]+1] + GNormalized[i,1:i-1].T
            p_xminusy[l[i:]+1] = p_xminusy[l[i:]+1] + GNormalized[i,i:s2].T

            # Contrast (CO)
            CO = CO + np.sum(((np.power(i+1 - o),2) * GNormalized[i,:]))

            # Inverse difference moment or homogenity (IDM)
            IDM = IDM + np.sum(GNormalized[i,o]/( 1 + np.power(i+1 - o,2)))

            # Dissimilarity (DIS)
            DIS = DIS + np.sum(np.abs(i+1 - o) * GNormalized[i,:])

            # Autocorrelation (ACORR)
            if typeflag.corr
                ACORR = ACORR + np.sum((i+1) * o * GNormalized[i,o])

            # Inverse difference (INV)
            INV = INV + np.sum(GNormalized[i,o] / ( 1 + np.abs(i+1 - o)))

            # Inverse difference normalized (INVN)
            INVN = INVN + np.sum(GNormalized[i,o] / ( 1 + (np.abs(i+1-o) / (s1))))

            # Inverse difference moment normalized (IDMN)
            IDMN = IDMN + np.sum(GNormalized[i,o] / ( 1 + np.power((i+1-o) / (s1),2)))


        #for i = 1:s1
        for i in range(s1):

            # Cluster shade (CLS)
            CLS = CLS + (np.power(i+1 + o - u_x - u_y,3)) * GNormalized[i,o].T

            # Cluster prominence (CLP)
            CLP = CLP + (np.power(i+1 + o - u_x - u_y,4)) * GNormalized[i,o].T

            sigma_x  = sigma_x  + np.sum((np.power((i+1) - u_x, 2)) * GNormalized[i,o])
            sigma_y  = sigma_y  + np.sum((np.power((o) - u_y, 2)) * GNormalized[i,o])


        # Correlation (CORR)
        if typeflag.corr
            CORR = (ACORR - (u_x*u_y)) / (sqrt(sigma_x)*np.sqrt(sigma_y))


        # Summed average (SA)
        #SA = (2:2*s1)*p_xplusy
        SA = range(2:2*s1+1) * p_xplusy

        # Summed entropy (SE)
        if typeflag.entropy
            SE = -np.sum(p_xplusy * np.log(p_xplusy))

        # Summed Variance (SV)
        SV = (np.power(range(2:2*s1+1) - SE, 2)) * p_xplusy

        # Difference varience (DV)
        DV = np.power(range(0:s1), 2) * p_xminusy

        # Difference entropy (DE)
        if typeflag.entropy
            DE = - np.sum(p_xminusy * np.log(p_xminusy))

            HXY = H

            for i in range(s1):
                HXY1 = HXY1 - GNormalized[i,:] * np.log(p_x[i] * p_y)
                HXY2 = HXY2 - np.sum(p_x[i] * p_y * np.log(p_x[i] * p_y))

            Hx = - np.sum(p_x * np.log(p_x))
            Hy = - np.sum(p_y * np.log(p_y))

            # Information measure of correlation 1 (IMC1)
            # Information measure of correlation 2 (IMC2)
            IMC1 = (HXY - HXY1) / ( np.max([Hx, Hy]))
            IMC2 = np.power(1 - np.exp(-2*(HXY2 - HXY))),0.5)

        # Maximum probability (MAXP)
        MAXP = np.max(GNormalized[:])

        if typeflag.texture || typeflag.global:
            Out[n,:] = [ACORR , CO, CORR, CLP, CLS ,  DIS, ASM, H, IDM, MAXP,
                SSV, SA, SV, SE, DV, DE, IMC1, IMC2, INV,INVN, IDMN]
        elif typeflag.corr:
            if typeflag.entropy
                Out[n,:] = [ACORR, CORR, H, SE, DE, IMC1, IMC2]
            else:
                Out[n,:] = [ACORR, CORR, IMC1, IMC2]
        else:
            Out[n,:] = [H, SE, DE, IMC1, IMC2]


    return reshape(Out, [1,np.shape(Out)[0]*np.shape(Out)[1]])
