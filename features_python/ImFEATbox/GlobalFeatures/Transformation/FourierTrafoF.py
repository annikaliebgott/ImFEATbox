import numpy as np
from scipy.misc import imrotate
from scipy.stats import moment

def FourierTrafoF(I, typeflag):
"""
    # Input:     - I: A 2D image
    #            - typeflag: Struct of logicals to permit extracting features
    #              based on desired characteristics:
    #                   + typeflag.global: all features
    #                   + typeflag.transform: all features
    #                   + typeflag.moments: only features based on moments
    #                   + typeflag.corr: only features based on correlation
    #              default: all features are being extracted
    #              For more information see README.txt
    #
    #
    # Output:    - Out: A (1x300) vector containing 300 metrics calculated from
    #                   the Fourier transform of an image
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

    if 'typeflag' not in globals():
        typeflag = dict()
        typeflag['global'] = True
        typeflag['transform'] = True
        typeflag['moments'] = True
        typeflag['corr'] = True


    # change image to a square matrix by means of zero padding
    length_x, length_y = np.shape(I)

    if (length_x < length_y):
        I = padarray(I, [( length_y - length_x ),0],'post')
    elif (length_x > length_y):
        I = padarray(I, [0, ( length_x - length_y )],'post')
    elif (length_x == length_y):
        I = double( I )

    #reservation for variables
    feat_corr = np.zeros(90)
    zero_sum = np.zeros(3)
    zero_sum1 = np.zeros(3)
    zero_sum2 = np.zeros(3)
    zero_sum3 = np.zeros(3)
    zero_sum4 = np.zeros(3)
    power_total = np.zeros(3)
    power_total1 = np.zeros(3)
    power_total2 = np.zeros(3)
    power_total3 = np.zeros(3)
    power_total4 = np.zeros(3)
    m_2_M = np.zeros(3)
    m_4_M = np.zeros(3)
    m_2_M1 = np.zeros(3)
    m_4_M1 = np.zeros(3)
    m_2_M2 = np.zeros(3)
    m_4_M2 = np.zeros(3)
    m_2_M3 = np.zeros(3)
    m_4_M3 = np.zeros(3)
    m_2_M4 = np.zeros(3)
    m_4_M4 = np.zeros(3)
    sd_M = np.zeros(3)
    sd_M1 = np.zeros(3)
    sd_M2 = np.zeros(3)
    sd_M3 = np.zeros(3)
    sd_M4 = np.zeros(3)
    sd_F = np.zeros(3)
    sd_F1 = np.zeros(3)
    sd_F2 = np.zeros(3)
    sd_F3 = np.zeros(3)
    sd_F4 = np.zeros(3)
    r_M = np.zeros(3)
    r_M1 = np.zeros(3)
    r_M2 = np.zeros(3)
    r_M3 = np.zeros(3)
    r_M4 = np.zeros(3)
    r_F = np.zeros(3)
    r_F1 = np.zeros(3)
    r_F2 = np.zeros(3)
    r_F3 = np.zeros(3)
    r_F4 = np.zeros(3)
    tr_M = np.zeros(3)
    tr_M1 = np.zeros(3)
    tr_M2 = np.zeros(3)
    tr_M3 = np.zeros(3)
    tr_M4 = np.zeros(3)
    tr_F = np.zeros(3)
    tr_F1 = np.zeros(3)
    tr_F2 = np.zeros(3)
    tr_F3 = np.zeros(3)
    tr_F4 = np.zeros(3)
    max_eig_M = np.zeros(3)
    max_eig_M1 = np.zeros(3)
    max_eig_M2 = np.zeros(3)
    max_eig_M3 = np.zeros(3)
    max_eig_M4 = np.zeros(3)
    max_eig_F = np.zeros(3)
    max_eig_F1 = np.zeros(3)
    max_eig_F2 = np.zeros(3)
    max_eig_F3 = np.zeros(3)
    max_eig_F4 = np.zeros(3)
    mean_eig_M = np.zeros(3)
    mean_eig_M1 = np.zeros(3)
    mean_eig_M2 = np.zeros(3)
    mean_eig_M3 = np.zeros(3)
    mean_eig_M4 = np.zeros(3)
    mean_eig_F = np.zeros(3)
    mean_eig_F1 = np.zeros(3)
    mean_eig_F2 = np.zeros(3)
    mean_eig_F3 = np.zeros(3)
    mean_eig_F4 = np.zeros(3)

    ## calculate transformation for different orientations

    for z in range(0, 3):

        if z == 1:
            # rotate image by 90 degree
            I = imrotate(I,90, interp='bilinear')   #'bicubic'
        elif z == 2:
            # rotate image by -90 degree
            I = imrotate(I,-180, interp='bilinear')

        # deterimie the Fourier - Transform
        F = np.fft.fft2(I)

        # create translation invariance
        M = np.abs(np.power(F,2))

        # break the image into sectors

        # get the center of the image
        half = np.max(np.shape(I))/2

        I1 = I[:half, :half)
        I2 = I[:half,half+1:np.max(np.shape(I)))
        I3 = I[half+1:np.max(np.shape(I)), :half)
        I4 = I[half+1:np.max(np.shape(I)), half+1:np.max(np.shape(I)))

        # Fourier Transform
        F1 = np.fft.fft2(I1)
        F2 = np.fft.fft2(I2)
        F3 = np.fft.fft2(I3)
        F4 = np.fft.fft2(I4)

        # create translation invariance
        M1 = np.abs(np.power(F1,2))
        M2 = np.abs(np.power(F2,2))
        M3 = np.abs(np.power(F3,2))
        M4 = np.abs(np.power(F4,2))

        ## extract features

        # determine the moments
        if (typeflag.moments || typeflag.global || typeflag.transform)
            m_2_M[z] = np.mean(moment(M,2))
            m_4_M[z] = np.mean(moment(M,4))
            m_2_M1[z] = np.mean(moment(M1,2))
            m_4_M1[z] = np.mean(moment(M1,4))
            m_2_M2[z] = np.mean(moment(M2,2))
            m_4_M2[z] = np.mean(moment(M2,4))
            m_2_M3[z] = np.mean(moment(M3,2))
            m_4_M3[z] = np.mean(moment(M3,4))
            m_2_M4[z] = np.mean(moment(M4,2))
            m_4_M4[z] = np.mean(moment(M4,4))
        end

        if (typeflag['global'] || typeflag['transform']):
            # standard derivation
            sd_M[z] = np.std(M)
            sd_M1[z] = np.std(M1)
            sd_M2[z] = np.std(M2)
            sd_M3[z] = np.std(M3)
            sd_M4[z] = np.std(M4)

            sd_F[z] = np.std(F)
            sd_F1[z] = np.std(F1)
            sd_F2[z] = np.std(F2)
            sd_F3[z] = np.std(F3)
            sd_F4[z] = np.std(F4)

            # rank
            r_M[z] = np.rank(M)
            r_M1[z] = np.rank(M1)
            r_M2[z] = np.rank(M2)
            r_M3[z] = np.rank(M3)
            r_M4[z] = np.rank(M4)

            r_F[z] = np.rank(F)
            r_F1[z] = np.rank(F1)
            r_F2[z] = np.rank(F2)
            r_F3[z] = np.rank(F3)
            r_F4[z] = np.rank(F4)

            # trace
            tr_M[z] = np.trace(M)
            tr_M1[z] = np.trace(M1)
            tr_M2[z] = np.trace(M2)
            tr_M3[z] = np.trace(M3)
            tr_M4[z] = np.trace(M4)

            tr_F[z] = np.trace(F)
            tr_F1[z] = np.trace(F1)
            tr_F2[z] = np.trace(F2)
            tr_F3[z] = np.trace(F3)
            tr_F4[z] = np.trace(F4)

            # largest/smallest eingenvalue
            eig_M = np.linalg.eig(M)[0]
            eig_M1 = np.linalg.eig(M1)[0]
            eig_M2 = np.linalg.eig(M2)[0]
            eig_M3 = np.linalg.eig(M3)[0]
            eig_M4 = np.linalg.eig(M4)[0]

            eig_F = np.linalg.eig(F)[0]
            eig_F1 = np.linalg.eig(F1)[0]
            eig_F2 = np.linalg.eig(F2)[0]
            eig_F3 = np.linalg.eig(F3)[0]
            eig_F4 = np.linalg.eig(F4)[0]

            max_eig_M[z] = np.max(eig_M)
            max_eig_M1[z] = np.max(eig_M1)
            max_eig_M2[z] = np.max(eig_M2)
            max_eig_M3[z] = np.max(eig_M3)
            max_eig_M4[z] = np.max(eig_M4)

            max_eig_F[z] = np.max(eig_F)
            max_eig_F1[z] = np.max(eig_F1)
            max_eig_F2[z] = np.max(eig_F2)
            max_eig_F3[z] = np.max(eig_F3)
            max_eig_F4[z] = np.max(eig_F4)

            # np.mean
            mean_eig_M[z] = np.mean(eig_M)
            mean_eig_M1[z] = np.mean(eig_M1)
            mean_eig_M2[z] = np.mean(eig_M2)
            mean_eig_M3[z] = np.mean(eig_M3)
            mean_eig_M4[z] = np.mean(eig_M4)

            mean_eig_F[z] = np.mean(eig_F)
            mean_eig_F1[z] = np.mean(eig_F1)
            mean_eig_F2[z] = np.mean(eig_F2)
            mean_eig_F3[z] = np.mean(eig_F3)
            mean_eig_F4[z] = np.mean(eig_F4)

        # correlation
        if (typeflag['corr'] or typeflag['transform'] or typeflag['global']):
            corr12 = np.corrcoef(M1,M2)
            corr13 = np.corrcoef(M1,M3)
            corr14 = np.corrcoef(M1,M4)
            corr32 = np.corrcoef(M3,M2)
            corr42 = np.corrcoef(M4,M2)
            corr34 = np.corrcoef(M4,M3)
            feat_corr((z-1)*30+1:(z-1)*30+30) = [np.std(corr12), np.mean(corr12),
                np.max(corr12(:)), np.min(corr12(:)), np.count_nonzero(corr12),
                np.std(corr13), np.mean(corr13), np.max(corr13(:)), np.min(corr13(:)), np.count_nonzero(corr13),
                np.std(corr14), np.mean(corr14), np.max(corr14(:)), np.min(corr14(:)), np.count_nonzero(corr14),
                np.std(corr32), np.mean(corr32), np.max(corr32(:)), np.min(corr32(:)), np.count_nonzero(corr32),
                np.std(corr42), np.mean(corr42), np.max(corr42(:)), np.min(corr42(:)), np.count_nonzero(corr42),
                np.std(corr34), np.mean(corr34), np.max(corr34(:)), np.min(corr34(:)), np.count_nonzero(corr34)]

        if (typeflag['global'] or typeflag['transform']):
            # count zero points
            zero_sum[z] = np.count_nonzero(F)
            zero_sum1[z] = np.count_nonzero(F1)
            zero_sum2[z] = np.count_nonzero(F2)
            zero_sum3[z] = np.count_nonzero(F3)
            zero_sum4[z] = np.count_nonzero(F4)


            # compute the DFT power

            # Transform length
            n = (2**ceil(log2(np.max(np.shape(I)))))**2

            # Power of the DFT
            power = np.power(F,conj(F))
            power1 = np.power(F1,conj(F1))
            power2 = np.power(F2,conj(F2))
            power3 = np.power(F3,conj(F3))
            power4 = np.power(F4,conj(F4))
            power_total[z] = np.sum( np.sum(power) ) /  (n*n)
            power_total1[z] = np.sum( np.sum(power1) ) / ( n*n )
            power_total2[z] = np.sum( np.sum(power2) ) / ( n*n )
            power_total3[z] = np.sum( np.sum(power3) ) / ( n*n )
            power_total4[z] = np.sum( np.sum(power4) ) / ( n*n )

    ## return feature vector
    if (typeflag['global'] or typeflag['transform']):
        return [zero_sum, zero_sum1, zero_sum2, zero_sum3, zero_sum4,
            power_total, power_total1, power_total2, power_total3, power_total4,
            m_2_M, m_4_M, m_2_M1, m_4_M1, m_2_M2, m_4_M2, m_2_M3, m_4_M3, m_2_M4, m_4_M4,
            sd_M, sd_M1, sd_M2, sd_M3, sd_M4,
            sd_F, sd_F1, sd_F2, sd_F3, sd_F4,
            r_M, r_M1, r_M2, r_M3, r_M4,
            r_F, r_F1, r_F2, r_F3, r_F4,
            tr_M, tr_M1, tr_M2, tr_M3, tr_M4,
            tr_F, tr_F1, tr_F2, tr_F3, tr_F4,
            max_eig_M, max_eig_M1, max_eig_M2, max_eig_M3, max_eig_M4,
            max_eig_F, max_eig_F1, max_eig_F2, max_eig_F3, max_eig_F4,
            mean_eig_M, mean_eig_M1, mean_eig_M2, mean_eig_M3, mean_eig_M4,
            mean_eig_F, mean_eig_F1, mean_eig_F2, mean_eig_F3, mean_eig_F4,
            feat_corr]
    elif not typeflag['moments']:
        return feat_corr
    elif not typeflag['corr']:
        return [m_2_M, m_4_M, m_2_M1, m_4_M1, m_2_M2, m_4_M2, m_2_M3, m_4_M3, m_2_M4, m_4_M4]
    else:
        return [m_2_M, m_4_M, m_2_M1, m_4_M1, m_2_M2, m_4_M2, m_2_M3, m_4_M3, m_2_M4, m_4_M4,
            feat_corr]
