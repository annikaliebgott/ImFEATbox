import numpy as np
from scipy.misc import imrotate
from ImFEATbox.__helperCommands import conv2float
from scipy.stats import skew, kurtosis

def SVDF(I, returnShape=False):
    """
     Input:     - I: A 2D image


     Output:    - Out: A (1x780) vector containing 780 metrics calculated
                  from singular value decomposition
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

    if returnShape:
        return (780,1)

    ## Calculate Singular Value Decomposition of the image

    # convert image to float
    I = conv2float(I)

    I = I.T

    # initialize feature variables
    dia_elements = np.zeros((np.min(np.shape(I)),3))
    eig_U = np.zeros((np.shape(I)[0],3))
    eig_V = np.zeros((np.shape(I)[1],3))
    det_U = np.zeros(3)
    det_V = np.zeros(3)
    trace_U = np.zeros(3)
    trace_V = np.zeros(3)
    rank_U = np.zeros(3)
    rank_V = np.zeros(3)
    median_eig_U = np.zeros(3)
    median_eig_V = np.zeros(3)
    max_eig_U = np.zeros(3)
    max_eig_V = np.zeros(3)
    mean_U = np.zeros(3)
    mean_V = np.zeros(3)
    mean_S = np.zeros(3)
    std_U = np.zeros(3)
    std_V = np.zeros(3)
    std_S = np.zeros(3)
    skewness_U = np.zeros(3)
    skewness_V = np.zeros(3)
    kurtosis_U = np.zeros(3)
    kurtosis_V = np.zeros(3)


    # Calculate the measures for 3 different orientations
    for z in range(0,3):

        if z == 1:
            # rotate image by 90 degree
            I = imrotate(I, 90, interp='bilinear')
        elif z == 2:
            # rotate image by -90 degree
            I = imrotate(I, -180, interp='bilinear')

        # calculate singular value decomposition with diagonal matrix S and
        # unitary matrices U and V
        [U,S,V] = np.linalg.svd(I)

        #U, V = U.T, V.T
        ## feature extraction

        # calculate diagonal elements of matrix S
        #for i in range(0, np.count_nonzero(S)):
        dia_elements[:,z] = S[:]

        # eigen values of U and V
        eig_U[:,z] = np.linalg.eig(U)[0]
        eig_V[:,z] = np.linalg.eig(V)[0]

        # determinant of U and V
        det_U[z] = np.linalg.det(U)
        det_V[z] = np.linalg.det(V)

        # trace of U and V
        trace_U[z] = np.trace(U)
        trace_V[z] = np.trace(V)

        # rank of U and V
        rank_U[z] = np.linalg.matrix_rank(U)
        rank_V[z] = np.linalg.matrix_rank(V)

        # skewness of U and V
        skewness_U[z] = skew(np.ndarray.flatten(U))
        skewness_V[z] = skew(np.ndarray.flatten(V))

        # kurtosis of U and V
        kurtosis_U[z] = kurtosis(np.ndarray.flatten(U), fisher=False, bias=False)
        kurtosis_V[z] = kurtosis(np.ndarray.flatten(V), fisher=False, bias=False)

        # mean of U, V and S
        mean_U[z] = np.mean(U)
        mean_V[z] = np.mean(V)
        mean_S[z] = np.mean(S)

        # standard deviation of U, V and S
        std_U[z] = np.std(U, ddof=1)
        std_V[z] = np.std(V, ddof=1)
        std_S[z] = np.std(S, ddof=1)

        # median of eigen values of U and V
        median_eig_U[z] = np.median(eig_U[:,z])
        median_eig_V[z] = np.median(eig_V[:,z])

        # maximum of eigen values of U and V
        max_eig_U[z] = np.max(eig_U[:,z])
        max_eig_V[z] = np.max(eig_V[:,z])

    ## return feature vector

    #np.prod(np.shape(eig_U[:100,:]))

    Out = np.hstack([np.ndarray.flatten(dia_elements[:40,:]),
        np.ndarray.flatten(eig_U[:100,:]),
        np.ndarray.flatten(eig_V[:100,:]),
        det_U, det_V, trace_U, trace_V, rank_U, rank_V, skewness_U, skewness_V,
        kurtosis_U, kurtosis_V, mean_U, mean_V, mean_S, std_U, std_V, std_S,
        median_eig_U, median_eig_V, max_eig_U, max_eig_V])

    return Out
