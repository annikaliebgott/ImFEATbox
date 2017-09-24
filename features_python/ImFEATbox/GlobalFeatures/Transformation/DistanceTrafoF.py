import numpy as np
from skimage.filters import threshold_otsu
from scipy import ndimage

def DistanceTrafoF(I,typeflag):
    """
     Input:     - I: A 2D image
                - typeflag: Struct of logicals to permit extracting features
                  based on desired characteristics:
                       + typeflag.global: all features
                       + typeflag.transform: all features
                       + typeflag.corr: only features based on correlation
                  default: all features are being extracted
                  For more information see README.txt


     Output:    - Out: A (1x56) vector containing 56 metrics calculated from
                  the distance transform
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

    #if ~exist('typeflag','var')
    if 'typeflag' not in globals():
        typeflag = dict()
        typeflag['global'] = True
        typeflag['transform'] = True
        typeflag['corr'] = True

    # converte image to real value float
    # im2bw can't process complex input values
    I = np.array(np.real(I), dtype='float')

    # threshold_otsu(): Global image threshold using Otsu's method
    BW = I > threshold_otsu(I)

    ## transfom image
    #[D1, IDX1] = bwdist(BW, 'chessboard')
    D1, IDX1 = ndimage.distance_transform_edt(BW, metric='chessboard', return_distances=True, return_indices=True)
    #[D2, IDX2] = bwdist(BW, 'cityblock') # taxicab = city block
    D2, IDX2 = ndimage.distance_transform_edt(BW, metric='taxicab', return_distances=True, return_indices=True)
    #[D3, IDX3] = bwdist(BW, 'euclidean')
    D3, IDX3 = ndimage.morphology.distance_transform_edt(BW, return_distances=True, return_indices=True) # euclidean

    # TODO find quasi-euclidean in python
    #[D4, IDX4] = bwdist(BW, 'quasi-euclidean')


    ## feature extraction

    # 2D correlation coefficient
    r12 = np.corrcoef(IDX1,IDX2)
    r13 = np.corrcoef(IDX1,IDX3)
    r14 = np.corrcoef(IDX1,IDX4)
    r23 = np.corrcoef(IDX2,IDX3)
    r24 = np.corrcoef(IDX2,IDX4)
    r34 = np.corrcoef(IDX3,IDX4)

    rD12 = np.corrcoef(D1,D2)
    rD13 = np.corrcoef(D1,D3)
    rD14 = np.corrcoef(D1,D4)
    rD23 = np.corrcoef(D2,D3)
    rD24 = np.corrcoef(D2,D4)
    rD34 = np.corrcoef(D3,D4)


    if (typeflag['global'] || typeflag['transform']):
        # mean of matrix elements
        B1 = np.mean(IDX1)
        B2 = np.mean(IDX2)
        B3 = np.mean(IDX3)
        B4 = np.mean(IDX4)

        B11 = np.mean(D1)
        B22 = np.mean(D2)
        B33 = np.mean(D3)
        B44 = np.mean(D4)

        # standard deviation of matrix elements
        S1 = np.std(IDX1, ddof=1)
        S2 = np.std(IDX2, ddof=1)
        S3 = np.std(IDX3, ddof=1)
        S4 = np.std(IDX4, ddof=1)

        S11 = np.std(D1, ddof=1)
        S22 = np.std(D2, ddof=1)
        S33 = np.std(D3, ddof=1)
        S44 = np.std(D4, ddof=1)

        s1 = np.std(np.std(IDX1, ddof=1), ddof=1)
        s2 = np.std(np.std(IDX2, ddof=1), ddof=1)
        s3 = np.std(np.std(IDX3, ddof=1), ddof=1)
        s4 = np.std(np.std(IDX4, ddof=1), ddof=1)

        s11 = np.std(np.std(D1, ddof=1), ddof=1)
        s22 = np.std(np.std(D2, ddof=1), ddof=1)
        s33 = np.std(np.std(D3, ddof=1), ddof=1)
        s44 = np.std(np.std(D4, ddof=1), ddof=1)

        # number of non zero elements
        nn_D11 = np.count_nonzero(IDX1)
        nn_D22 = np.count_nonzero(IDX2)
        nn_D33 = np.count_nonzero(IDX3)
        nn_D44 = np.count_nonzero(IDX4)

        nn_D1 = np.count_nonzero(D1)
        nn_D2 = np.count_nonzero(D2)
        nn_D3 = np.count_nonzero(D3)
        nn_D4 = np.count_nonzero(D4)

        # determine max/min values of IDX
        max_IDX1 = np.max(IDX1)
        max_IDX2 = np.max(IDX2)
        max_IDX3 = np.max(IDX3)
        max_IDX4 = np.max(IDX4)
        min_IDX1 = np.min(IDX1)
        min_IDX2 = np.min(IDX2)
        min_IDX3 = np.min(IDX3)
        min_IDX4 = np.min(IDX4)

        # determine max value of D
        max_D1 = np.max(D1)
        max_D2 = np.max(D2)
        max_D3 = np.max(D3)
        max_D4 = np.max(D4)

    ## return feature vector

    if not (typeflag['global'] || typeflag['transform']):
        return [r12, r13, r14, r23, r24, r34, rD12, rD13, rD14, rD23, rD24, rD34]
    else:
        return [B1, B2, B3, B4, B11, B22, B33, B44,
            S1, S2, S3, S4, S11, S22, S33, S44,
            s1, s2, s3, s4, s11, s22, s33, s44,
            r12, r13, r14, r23, r24, r34, rD12, rD13, rD14, rD23, rD24, rD34,
            nn_D1, nn_D2, nn_D3, nn_D4, nn_D11, nn_D22, nn_D33, nn_D44,
            max_IDX1, max_IDX2, max_IDX3, max_IDX4, min_IDX1, min_IDX2, min_IDX3, min_IDX4,
            max_D1, max_D2, max_D3, max_D4]
