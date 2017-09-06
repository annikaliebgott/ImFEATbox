import numpy as np
from ImFEATbox.__helperCommands import rgb2grayscale, isColorImage, conv2float
from scipy.signal import convolve2d

def GradientF(I, typeflag=None, gradtype=None, returnShape=False):
    """
         Input:     - I: A 2D image
                    - typeflag: Struct of logicals to permit extracting features
                      based on desired characteristics:
                           + typeflag['global']: all features
                           + typeflag['texture']: all features
                           + typeflag['gradient']: all features
                           + typeflag['entropy']: only features based on entropy
                      default: all features are being extracted
                      For more information see README.txt
                    - gradtype: Struct of logicals to choose which type of
                      gradient should be applied:
                           + gradtype['first']: first order gradient
                           + gradtype['second']: second order gradient (Laplacian)
                      default: both types of gradients are used
                    - returnShape: A logical flag to return only the shape of features.
                      Default: False


         Output:    - Out: A (1x81) vector containing 81 metrics calculated from
                      image gradients
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
    # Contact: annika.liebgott@iss.uni-stuttgart.de
    # ************************************************************************
    #
    # Implementation based on:  McGee et al.: "Image metric-based correction
    #                           (Autocorrection) of motion effects: Analysis of
    #                           image metrics", Journal of Magnetic Resonance
    #                           Imaging, vol. 11, issue 2, p. 174-181, Feb 2000

    ##




    if typeflag == None:
        typeflag = dict()
        typeflag['global'] = True
        typeflag['texture'] = True
        typeflag['entropy'] = True


    if gradtype == None:
        gradtype = dict()
        gradtype['first'] = True
        gradtype['second'] = True



    if returnShape:
        ## returns only the shape of feature vectors depending on parameters
        if typeflag['global'] or typeflag['texture'] or typeflag['gradient']:
            if gradtype['first'] and gradtype['second']:
                return (81,1)
            elif gradtype['first']:
                return (36,1)
            else:
                return (45,1)
        elif gradtype['first'] and gradtype['second']:
            return (9,1)
        elif gradtype['first']:
            return (4,1)
        else:
            return (5,1)


    # Check for color image and convert to grayscale
    if isColorImage(I):
        I = rgb2grayscale(I)

    # make float
    I = conv2float(I)

    Height, Width = np.shape(I)
    n = Height * Width

    ## define gradients
    if gradtype['first']:
        # gradient 1
        # gradient direction: x
        g1_x = np.array([   [1, -1],
                            [0, 0]])
        # gradient direction: y
        g1_y = np.array([   [1, 0],
                            [-1, 0]])

        # gradient 2
        g2_x = np.array([   [0, 0, 0],
                            [1, 0, -1],
                            [0, 0, 0]])

        g2_y = np.array([   [0, 1, 0],
                            [0, 0, 0],
                            [0, -1, 0]])

    if gradtype['second']:
        # laplacian 1
        l1_x = np.array([   [0, 0, 0],
                            [-1, 2, -1],
                            [0, 0, 0]])

        l1_y = np.array([   [0, -1, 0],
                            [0, 2, 0],
                            [0, -1, 0]])

        # laplacian 2
        l2 = np.array([ [-1, -2, -1],
                        [-2, 12, -2],
                        [-1, -2, -1]])

        # laplacian 3
        l3 = np.array([ [0, -1, 0],
                        [-1, 4, -1],
                        [0, -1, 0]])

        # laplacian 4
        l4 = np.array([ [-1, -1, -1],
                        [-1, 8, -1],
                        [-1, -1, -1]])

    ## convolve image
    if gradtype['first']:
        G1_x = convolve2d(I,g1_x, boundary='fill', mode='full')
        G1_y = convolve2d(I,g1_y, boundary='symm', mode='full')
        G2_x = convolve2d(I,g2_x, boundary='symm', mode='full')
        G2_y = convolve2d(I,g2_y, boundary='symm', mode='full')

    if gradtype['second']:
        L1_x = convolve2d(I,l1_x, boundary='symm', mode='full')
        L1_y = convolve2d(I,l1_y, boundary='symm', mode='full')
        L2 = convolve2d(I,l2, boundary='symm', mode='full')
        L3 = convolve2d(I,l3, boundary='symm', mode='full')
        L4 = convolve2d(I,l4, boundary='symm', mode='full')

    ## extract features

    # summed gradients
    if gradtype['first']:
        G1_y_sum = np.sum(np.sum(np.abs(G1_y)))
        G1_x_sum = np.sum(np.sum(np.abs(G1_x)))
        G2_x_sum = np.sum(np.sum(np.abs(G2_x)))
        G2_y_sum = np.sum(np.sum(np.abs(G2_y)))

    if gradtype['second']:
        L1_x_sum = np.sum(np.sum(np.abs(L1_x)))
        L1_y_sum = np.sum(np.sum(np.abs(L1_y)))
        L2_sum = np.sum(np.sum(np.abs(L2)))
        L3_sum = np.sum(np.sum(np.abs(L3)))
        L4_sum = np.sum(np.sum(np.abs(L4)))

    # normalized gradients
    if gradtype['first']:
        G1_x_norm = np.abs(G1_x)/G1_x_sum
        G1_y_norm = np.abs(G1_y)/G1_y_sum
        G2_x_norm = np.abs(G2_x)/G2_x_sum
        G2_y_norm = np.abs(G2_y)/G2_y_sum

    if gradtype['second']:
        L1_x_norm = np.abs(L1_x)/L1_x_sum
        L1_y_norm = np.abs(L1_y)/L1_y_sum
        L2_norm = np.abs(L2)/L2_sum
        L3_norm = np.abs(L3)/L3_sum
        L4_norm = np.abs(L4)/L4_sum

    if typeflag['global'] or typeflag['texture'] or typeflag['gradient']:
        # sum of squared gradients
        if gradtype['first']:
            G1_x_2 = np.sum(np.sum(np.power(G1_x,2)))/n
            G1_y_2 = np.sum(np.sum(np.power(G1_y,2)))/n
            G2_x_2 = np.sum(np.sum(np.power(G2_x,2)))/n
            G2_y_2 = np.sum(np.sum(np.power(G2_y,2)))/n
        if gradtype['second']:
            L1_x_2 = np.sum(np.sum(np.power(L1_x,2)))/n
            L1_y_2 = np.sum(np.sum((L1_y,2)))/n
            L2_2 = np.sum(np.sum((L2,2)))/n
            L3_2 = np.sum(np.sum((L3,2)))/n
            L4_2 = np.sum(np.sum((L4,2)))/n

        # sum of 4th power of gradients
        if gradtype['first']:
            G1_x_4 = np.sum(np.sum(np.power(G1_x,4)))/n
            G1_y_4 = np.sum(np.sum(np.power(G1_y,4)))/n
            G2_x_4 = np.sum(np.sum(np.power(G2_x,4)))/n
            G2_y_4 = np.sum(np.sum(np.power(G2_y,4)))/n

        if gradtype['second']:
            L1_x_4 = np.sum(np.sum(np.power(L1_x,4)))/n
            L1_y_4 = np.sum(np.sum(np.power(L1_y,4)))/n
            L2_4 = np.sum(np.sum(np.power(L2,4)))/n
            L3_4 = np.sum(np.sum(np.power(L3,4)))/n
            L4_4 = np.sum(np.sum(np.power(L4,4)))/n

        # maximum of normalized gradients
        if gradtype['first']:
            G1_x_norm_max = np.max(G1_x_norm[:])
            G1_y_norm_max = np.max(G1_y_norm[:])
            G2_x_norm_max = np.max(G2_x_norm[:])
            G2_y_norm_max = np.max(G2_y_norm[:])

        if gradtype['second']:
            L1_x_norm_max = np.max(L1_x_norm[:])
            L1_y_norm_max = np.max(L1_y_norm[:])
            L2_norm_max = np.max(L2_norm[:])
            L3_norm_max = np.max(L3_norm[:])
            L4_norm_max = np.max(L4_norm[:])

        # standard deviation of normalized gradients
        if gradtype['first']:
            G1_x_norm_std = np.std(G1_x_norm)
            G1_y_norm_std = np.std(G1_y_norm)
            G2_x_norm_std = np.std(G2_x_norm)
            G2_y_norm_std = np.std(G2_y_norm)

        if gradtype['second']:
            L1_x_norm_std = np.std(L1_x_norm)
            L1_y_norm_std = np.std(L1_y_norm)
            L2_norm_std = np.std(L2_norm)
            L3_norm_std = np.std(L3_norm)
            L4_norm_std = np.std(L3_norm)

        # mean of normalized gradients
        if gradtype['first']:
            G1_x_norm_mean = np.mean(G1_x_norm)
            G1_y_norm_mean = np.mean(G1_y_norm)
            G2_x_norm_mean = np.mean(G2_x_norm)
            G2_y_norm_mean = np.mean(G2_y_norm)

        if gradtype['second']:
            L1_x_norm_mean = np.mean(L1_x_norm)
            L1_y_norm_mean = np.mean(L1_y_norm)
            L2_norm_mean = np.mean(L2_norm)
            L3_norm_mean = np.mean(L3_norm)
            L4_norm_mean = np.mean(L3_norm)

        # sum of normalized gradients squared
        if gradtype['first']:
            G1_x_norm_2 = np.sum(np.sum(np.power(G1_x_norm,2)))/n
            G1_y_norm_2 = np.sum(np.sum(np.power(G1_y_norm,2)))/n
            G2_x_norm_2 = np.sum(np.sum(np.power(G2_x_norm,2)))/n
            G2_y_norm_2 = np.sum(np.sum(np.power(G2_y_norm,2)))/n

        if gradtype['second']:
            L1_x_norm_2 = np.sum(np.sum(np.power(L1_x_norm,2)))/n
            L1_y_norm_2 = np.sum(np.sum(np.power(L1_y_norm,2)))/n
            L2_norm_2 = np.sum(np.sum(np.power(L2_norm,2)))/n
            L3_norm_2 = np.sum(np.sum(np.power(L3_norm,2)))/n
            L4_norm_2 = np.sum(np.sum(np.power(L4_norm,2)))/n

        # sum of normalized gradients to 4th power
        if gradtype['first']:
            G1_x_norm_4 = np.sum(np.sum(np.power(G1_x_norm,4)))/n
            G1_y_norm_4 = np.sum(np.sum(np.power(G1_y_norm,4)))/n
            G2_x_norm_4 = np.sum(np.sum(np.power(G2_x_norm,4)))/n
            G2_y_norm_4 = np.sum(np.sum(np.power(G2_y_norm,4)))/n

        if gradtype['second']:
            L1_x_norm_4 = np.sum(np.sum(np.power(L1_x_norm,4)))/n
            L1_y_norm_4 = np.sum(np.sum(np.power(L1_y_norm,4)))/n
            L2_norm_4 = np.sum(np.sum(np.power(L2_norm,4)))/n
            L3_norm_4 = np.sum(np.sum(np.power(L3_norm,4)))/n
            L4_norm_4 = np.sum(np.sum(np.power(L4_norm,4)))/n

    # marginal entropies of normalized gradients
    if gradtype['first']:
        G1_x_E = -np.sum(G1_x_norm[np.nonzero(G1_x_norm)] * np.log2(G1_x_norm[np.nonzero(G1_x_norm)]))
        G1_y_E = -np.sum(G1_y_norm[G1_y_norm != 0] * np.log2(G1_y_norm[G1_y_norm != 0]))
        G2_x_E = -np.sum(G2_x_norm[G2_x_norm != 0] * np.log2(G2_x_norm[G2_x_norm != 0]))
        G2_y_E = -np.sum(G2_y_norm[G2_y_norm != 0] * np.log2(G2_y_norm[G2_y_norm != 0]))

    if gradtype['second']:
        L1_x_E = -np.sum(L1_x_norm[L1_x_norm != 0] * np.log2(L1_x_norm[L1_x_norm != 0]))
        L1_y_E = -np.sum(L1_y_norm[L1_y_norm != 0] * np.log2(L1_y_norm[L1_y_norm != 0]))
        L2_E = -np.sum(L2_norm[L2_norm != 0] * np.log2(L2_norm[L2_norm != 0]))
        L3_E = -np.sum(L3_norm[L3_norm != 0] * np.log2(L3_norm[L3_norm != 0]))
        L4_E = -np.sum(L4_norm[L4_norm != 0] * np.log2(L4_norm[L4_norm != 0]))

    ## Build output vector Out
    if typeflag['global'] or typeflag['texture'] or typeflag['gradient']:
        if gradtype['first'] and gradtype['second']:
            Out = np.array([G1_x_sum,G1_y_sum,G2_x_sum,G2_y_sum,L1_x_sum,L1_y_sum,L2_sum,L3_sum,L4_sum,
                G1_x_2,G1_y_2,G2_x_2,G2_y_2,L1_x_2,L1_y_2,L2_2,L3_2,L4_2,
                G1_x_4,G1_y_4,G2_x_4,G2_y_4,L1_x_4,L1_y_4,L2_4,L3_4,L4_4,
                G1_x_norm_max,G1_y_norm_max,G2_x_norm_max,G2_y_norm_max,
                L1_x_norm_max,L1_y_norm_max,L2_norm_max,L3_norm_max,L4_norm_max,
                G1_x_norm_std,G1_y_norm_std,G2_x_norm_std,G2_y_norm_std,
                L1_x_norm_std,L1_y_norm_std,L2_norm_std,L3_norm_std,L4_norm_std,
                G1_x_norm_mean,G1_y_norm_mean,G2_x_norm_mean,G2_y_norm_mean,
                L1_x_norm_mean,L1_y_norm_mean,L2_norm_mean,L3_norm_mean,L4_norm_mean,
                G1_x_norm_2,G1_y_norm_2,G2_x_norm_2,G2_y_norm_2,
                L1_x_norm_2,L1_y_norm_2,L2_norm_2,L3_norm_2,L4_norm_2,
                G1_x_norm_4,G1_y_norm_4,G2_x_norm_4,G2_y_norm_4,
                L1_x_norm_4,L1_y_norm_4,L2_norm_4,L3_norm_4,L4_norm_4,
                G1_x_E,G1_y_E,G2_x_E,G2_y_E,L1_x_E,L1_y_E,L2_E,L3_E,L4_E])
        elif gradtype['first']:
            Out = [G1_x_sum,G1_y_sum,G2_x_sum,G2_y_sum,
                G1_x_2,G1_y_2,G2_x_2,G2_y_2,
                G1_x_4,G1_y_4,G2_x_4,G2_y_4,
                G1_x_norm_max,G1_y_norm_max,G2_x_norm_max,G2_y_norm_max,
                G1_x_norm_std,G1_y_norm_std,G2_x_norm_std,G2_y_norm_std,
                G1_x_norm_mean,G1_y_norm_mean,G2_x_norm_mean,G2_y_norm_mean,
                G1_x_norm_2,G1_y_norm_2,G2_x_norm_2,G2_y_norm_2,
                G1_x_norm_4,G1_y_norm_4,G2_x_norm_4,G2_y_norm_4,
                G1_x_E,G1_y_E,G2_x_E,G2_y_E]
        else:
            Out = [L1_x_sum,L1_y_sum,L2_sum,L3_sum,L4_sum,
                L1_x_2,L1_y_2,L2_2,L3_2,L4_2,
                L1_x_4,L1_y_4,L2_4,L3_4,L4_4,
                L1_x_norm_max,L1_y_norm_max,L2_norm_max,L3_norm_max,L4_norm_max,
                L1_x_norm_std,L1_y_norm_std,L2_norm_std,L3_norm_std,L4_norm_std,
                L1_x_norm_mean,L1_y_norm_mean,L2_norm_mean,L3_norm_mean,L4_norm_mean,
                L1_x_norm_2,L1_y_norm_2,L2_norm_2,L3_norm_2,L4_norm_2,
                L1_x_norm_4,L1_y_norm_4,L2_norm_4,L3_norm_4,L4_norm_4,
                L1_x_E,L1_y_E,L2_E,L3_E,L4_E]
    elif gradtype['first'] and gradtype['second']:
        Out = [G1_x_E,G1_y_E,G2_x_E,G2_y_E,L1_x_E,L1_y_E,L2_E,L3_E,L4_E]
    elif gradtype['first']:
        Out = [G1_x_E,G1_y_E,G2_x_E,G2_y_E]
    else:
        Out = [L1_x_E,L1_y_E,L2_E,L3_E,L4_E]

    return Out
