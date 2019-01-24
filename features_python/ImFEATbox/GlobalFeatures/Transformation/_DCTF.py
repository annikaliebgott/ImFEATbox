import numpy as np
from scipy.misc import imrotate
from scipy.fftpack import dct
from numpy.linalg import eig, svd, det
from ImFEATbox.__helperCommands import conv2float

def DCTF(I, typeflag=None, returnShape=False):
    """
     Input:     - I: A 2D image
                - typeflag: Struct of logicals to permit extracting features
                  based on desired characteristics:
                       + typeflag['global']: all features
                       + typeflag['transform']: all features
                       + typeflag['corr']: only features based on correlation
                  default: all features are being extracted
                  For more information see README.txt

     Output:    - Out: A (1x2901) vector containing 2901 metrics calculated
                       from the discrete cosine transform
   """
    # ************************************************************************
    # Implemented for MRI feature extraction by the Department of Diagnostic
    # and Interventional Radiology, University Hospital of Tuebingen, Germany
    # and the Institute of Signal Processing and System Theory University of
    # Stuttgart, Germany. Last modified: December 2016
    #
    # This implementation is part of ImFEATbox, a toolbox for image feature
    # extraction and analysis. Available online at:
    # https://github.com/annikaliebgott/ImFEATbox
    #
    # Contact: annika.liebgott@iss.uni-stuttgart.de
    # ************************************************************************

    if returnShape:
        return (2901,1)

    if typeflag == None:
        typeflag = dict()
        typeflag['global'] = True
        typeflag['transform'] = True
        typeflag['corr'] = True


    I = conv2float(I)

    # reserve space for the feature vector
    # if typeflag['transform'] || typeflag['global'] == true: extract all features
    # if only typeflag['corr'] == true: extract especially correlation based
    # features
    if typeflag['transform'] or typeflag['global']:
        f = np.zeros((3,967))
    else:
        f = np.zeros((3,84))

    ## calculate DCT for different orientations: horizontal,vertical, diagonal

    for z in range(0, 3):
        # horizontal, vertical and diagonal by counting the number of occurrence
        # of each selected coefficient with respect to its position
        # equal importance of horizontal, vertical and diagonal

        if (z == 1):
            # rotate image by 90 degree
            I = imrotate(I, 90, interp='bilinear')  #'bicubic'
        elif (z == 2):
            # rotate image by -90 degree
            I = imrotate(I, -180, interp='bilinear')

        # index counter
        idx = 0

        # 2D discrete cosine transform, later on used as reference matrix
        B = dct(I, type=2, norm='ortho')

        # perform SVD
        U,S,V = svd(B)
        eUB = eig(U)[0]
        eVB = eig(V)[0]

        # extract features
        if typeflag['transform'] or typeflag['global']:
            f[z,idx:idx+99] = [np.swapaxes(eUB[:50], 0, 1), np.swapaxes(eVB[:50], 0, 1)]
            idx = idx+99

            if np.shape(B)[1] < 150:
                coefB = [B[0,0], B[0,39], B[19,59], B[79,99], B[99, np.shape(B)[1]]]
            else:
                coefB = [B[0,0], B[0,39], B[19,59], B[79,99], B[99, 149]]

            f[z,idx:idx + len(coefB)-1] = coefB
            idx = idx + len(coefB)

            f[z,idx:idx+3] = [det(U), trace(U), det(V), trace(V)]
            idx = idx + 4


        ## calculate transform for different block decomposition sizes
        for i in [2, 4, 8, 16, 32, 64]:

            # TODO block_struct ?!

            # extract spatio-temporal features
            # decompose blocks of images into independent tiles
            #fun = @(block_struct) np.std(block_struct.data, ddof=1) * np.ones(np.shape(block_struct.data))
            #I2 = blockproc(I,[i i],fun)
            I2 = 0

            # 2D discrete cosine transform
            # I2 has to be grayscale
            B2 = dct(I2, type=2, norm='ortho')

            # output of some coefficients
            # 3D zigzag transversal to select 25# of the coefficents
            if typeflag['transform'] or typeflag['global']:
                if np.shape(B2)[1] < 150:
                    coefB2 = [B[0,0], B[0,39], B[19,59], B[79,99], B[99, np.shape(B)[1]]]
                else:
                    coefB2 = [B[0,0], B[0,39], B[19,59], B[79,99], B[99, 149]]
                end
                f[z,idx+1:idx + len(coefB2)] = coefB2
                idx = idx + len(coefB2)

                # calculate important features of the 2D DCT
                temp = [np.std(B2, ddof=1), np.std(np.std(B2, ddof=1), ddof=1), np.mean(B2), np.linalg.matrix_rank(B2), np.max(B2), np.min(B2), np.count_nonzero(B2)]
                f[z,idx+1:idx+ len(temp)] = temp
                idx = idx + len(temp)
            end

            # correlation between DCT of I and decomposed I
            Corr = np.corrcoef(B, B2)
            f[z,idx] = Corr(1,2)
            idx = idx +1

            # singular value deomposition, U & V are square matrices
            [U2,S2,V2] = svd(B2)
            eUB2 = eig(U2)[0]
            eVB2 = eig(V2)[0]

            if typeflag['transform'] or typeflag['global']:
                temp = [np.swapaxes(eUB2[:50], 0, 1), np.swapaxes(eVB2[:50], 0, 1),
                np.std(U2, ddof=1), np.std(S2, ddof=1), np.std(V2, ddof=1), np.mean(U2), np.mean(S2), np.mean(V2),
                np.max(U2), np.max(S2), np.max(V2), np.min(U2), np.min(S2), np.min(V2),
                np.count_nonzero(B2), det(U2), det(V2), trace(U2), trace(V2)]
                f[z,idx:idx + len(temp)-1] = temp
                idx = idx + len(temp)

            # correlation between SVD matrices
            corrS = np.corrcoef(S,S2)
            corrU = np.corrcoef(U,U2)
            corrV = np.corrcoef(V,V2)

            # correlation between eigenvalues of SVD
            correV = np.corrcoef(eVB,eVB2)
            correU = np.corrcoef(eUB,eUB2)
            xcorreV = np.correlate(eVB,eVB2)
            xcorreU = np.correlate(eUB,eUB2)
            temp = [corrS[0,1], corrU[0,1], corrV[0,1], correV[0,1], correU[0,1], mean(xcorreV), std(xcorreV, ddof=1), np.min(xcorreV), np.max(xcorreV). np.mean(xcorreU), std(xcorreU, ddof=1), np.min(xcorreU), np.max(xcorreU)]

            f[z,idx:idx + len(temp) -1] = temp
            idx = idx + len(temp)
        end
    end

    ## return feature vector
    f = np.swapaxes(f, 0, 1)
    Out = f #.T

    return Out
