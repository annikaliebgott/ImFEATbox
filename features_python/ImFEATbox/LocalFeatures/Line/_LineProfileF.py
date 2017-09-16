import numpy as np
from scipy.stats import moment
from skimage.measure import profile_line

def LineProfileF(I, typeflag=None, plotflag=False, returnShape=False):
    """
 Input:     - I: A 2D image
            - typeflag: Struct of logicals to permit extracting features
              based on desired characteristics:
                   + typeflag['local']: all features
                   + typeflag['texture']: all features
                   + typeflag['corr']: only features based on correlation
                   + typeflag['moments']: only features based on moments
              default: all features are being extracted
              For more information see README.txt
            - plotflag: logical flag to enable/disable visualization.
              default: plotflag = false


 Output:    - Out: A (1x122) vector containing 122 metrics calculated from
                   an intensity profile of the image
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
        return (122,1)

    if typeflag == None:
        typeflag = dict()
        typeflag['local'] = True
        typeflag['texture'] = True
        typeflag['corr'] = True
        typeflag['moments'] = True

    if plotflag == None:
        plotflag = False

    ## Extract the intensity profile along different line segments

    x_max, y_max = np.shape(I)
    x_half = x_max / 2.0
    y_half = y_max / 2.0

    # diagonal lines with fitting size for the input image
    x1 = [0, x_max]
    y1 = [0, y_max]
    x2 = [x_max, 0]
    y2 = [0, y_max]

    # vertical line through the center of the image
    x3 = [x_half, x_half]
    y3 = [0, y_max]

    # horizontal line through the center of the image
    x4 = [0, x_max]
    y4 = [y_half, y_half]

    # view line profile
    if plotflag:
        x = [x1, x2]
        y = [y1, y2]
        profile_line(I, src=(0,0), dst=(x1[1],y2[1]))
        # TODO draw with matplotlib


    # extract intensity profiles along the defined lines
    c1 = profile_line(I, src=(x1[0],y1[0]), dst=(x1[1],y1[1]))
    c2 = profile_line(I, src=(x2[0],y2[0]), dst=(x2[1],y2[1]))
    c3 = profile_line(I, src=(x3[0],y3[0]), dst=(x3[1],y3[1]))
    c4 = profile_line(I, src=(x4[0],y4[0]), dst=(x4[1],y4[1]))

    l1 = np.shape(c1)[0]
    l2 = np.shape(c2)[0]
    l3 = np.shape(c3)[0]
    l4 = np.shape(c4)[0]

    # Use zero-padding to get the same size for all arrays
    m = np.max([l1, l2, l3, l4])
    c1 = np.pad(c1, (0,max(0,(m-l1)-len(c1))) ,mode='constant')
    c1 = np.pad(c2, (0,max(0,(m-l2)-len(c2))) ,mode='constant')
    c1 = np.pad(c3, (0,max(0,(m-l3)-len(c3))) ,mode='constant')
    c1 = np.pad(c4, (0,max(0,(m-l4)-len(c4))) ,mode='constant')

    ## Feature extraction

    # linear correlation
    if (typeflag['corr'] or typeflag['local'] or typeflag['texture']):
        print(np.shape(c1))
        print(np.shape(c2))
        roh12 = np.corrcoef(c1,c2)
        roh13 = np.corrcoef(c1,c3)
        roh14 = np.corrcoef(c1,c4)
        roh23 = np.corrcoef(c2,c3)
        roh24 = np.corrcoef(c4,c2)
        roh34 = np.corrcoef(c3,c4)

    # central moments
    if (typeflag['moments'] or typeflag['local'] or typeflag['texture']):
        m1 = moment(c1,2)
        m2 = moment(c2,2)
        m3 = moment(c3,2)
        m4 = moment(c4,2)

        m1_5 = moment(c1,5)
        m2_5 = moment(c2,5)
        m3_5 = moment(c3,5)
        m4_5 = moment(c4,5)

    if (typeflag['local'] or typeflag['texture']):
        # average value of array
        mean1 = np.mean(c1)
        mean2 = np.mean(c2)
        mean3 = np.mean(c3)
        mean4 = np.mean(c4)

        # standard deviation
        sd1 = np.std(c1)
        sd2 = np.std(c2)
        sd3 = np.std(c3)
        sd4 = np.std(c4)

        # Percentiles of a data set
        pr1 = np.percentile(c1,48)
        pr2 = np.percentile(c2,48)
        pr3 = np.percentile(c3,48)
        pr4 = np.percentile(c4,48)

        # fast fourier transform
        fft1 = np.fft.fft(c1)
        fft2 = np.fft.fft(c2)
        fft3 = np.fft.fft(c3)
        fft4 = np.fft.fft(c4)

        # power of spectrum
        p1 = np.abs(np.sum(np.power(fft1,2)) / len(fft1))
        p2 = np.abs(np.sum(np.power(fft2,2)) / len(fft2))
        p3 = np.abs(np.sum(np.power(fft3,2)) / len(fft3))
        p4 = np.abs(np.sum(np.power(fft4,2)) / len(fft4))

        # derivates
        dc1 = np.diff(c1)
        dc2 = np.diff(c2)
        dc3 = np.diff(c3)
        dc4 = np.diff(c4)

        # average value of array of the derivates
        mean_dc1 = np.mean(dc1)
        mean_dc2 = np.mean(dc2)
        mean_dc3 = np.mean(dc3)
        mean_dc4 = np.mean(dc4)

        # standard deviation of the derivates
        sd_dc1 = np.std(dc1)
        sd_dc2 = np.std(dc2)
        sd_dc3 = np.std(dc3)
        sd_dc4 = np.std(dc4)

        # fast fourier transform of the derivates
        fft_dc1 = np.fft.fft(dc1)
        fft_dc2 = np.fft.fft(dc2)
        fft_dc3 = np.fft.fft(dc3)
        fft_dc4 = np.fft.fft(dc4)

        # power of spectrum of the derivates
        p_dc1 = np.abs(np.sum(np.power(fft_dc1,2)) / len(fft_dc1))
        p_dc2 = np.abs(np.sum(np.power(fft_dc2,2)) / len(fft_dc2))
        p_dc3 = np.abs(np.sum(np.power(fft_dc3,2)) / len(fft_dc3))
        p_dc4 = np.abs(np.sum(np.power(fft_dc4,2)) / len(fft_dc4))


    ## Return feature vector

    if not typeflag['local'] or typeflag['texture']:
        if not typeflag['moments']:
            Out = np.array([roh12, roh14, roh24, roh13, roh23, roh34])
        elif not typeflag['corr']:
            Out = np.array([m1, m2, m3, m4, m1_5, m2_5, m3_5, m4_5])
        else:
            Out = np.array([roh12, roh14, roh24, roh13, roh23, roh34,
                m1, m2, m3, m4, m1_5, m2_5, m3_5, m4_5])
    else:
        Out = np.array([roh12, roh14, roh24, roh13, roh23, roh34,
            mean1, mean2, mean3, mean4,
            sd1, sd2, sd3, sd4,
            pr1, pr2, pr3, pr4,
            m1, m2, m3, m4, m1_5, m2_5, m3_5, m4_5,
            p1, p2, p3, p4,
            fft1[:10,1], fft2[:10,1], fft3[:10,1], fft4[:10,1],
            mean_dc1, mean_dc2, mean_dc3, mean_dc4,
            sd_dc1, sd_dc2, sd_dc3, sd_dc4,
            p_dc2, p_dc1, p_dc3, p_dc4,
            fft_dc1[:10,1], fft_dc2[:10,1], fft_dc3[:10,1], fft_dc4[:10,1]])
    return Out
