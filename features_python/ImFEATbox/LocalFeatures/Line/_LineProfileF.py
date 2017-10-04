import numpy as np
import warnings
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


    if typeflag == None:
        typeflag = dict()
        typeflag['local'] = True
        typeflag['texture'] = True
        typeflag['corr'] = True
        typeflag['moments'] = True

    if typeflag['local']  == typeflag['texture'] == typeflag['corr'] == typeflag['moments'] == False:
        # catching this case. lineprofile like this does not make sense.
        # so we throw a warning and set all to True
        typeflag['local'] = True
        typeflag['texture'] = True
        typeflag['corr'] = True
        typeflag['moments'] = True
        warnings.warn("typeflag local, texture, corr and moments are false. Using default settings.")

    if returnShape:
        if not (typeflag['local'] or typeflag['texture']):
            if not typeflag['moments']:
                return (6,1)
            elif not typeflag['corr']:
                return (8,1)
            else:
                return (14,1)
        else:
            return (122,1)

    if plotflag == None:
        plotflag = False

    ## Extract the intensity profile along different line segments

    y_max, x_max = np.shape(I)
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
        profile_line(I, src=(0,0), dst=(y1[1],x2[1]))
        # TODO draw with matplotlib

    # extract intensity profiles along the defined lines
    c1 = profile_line(I, src=(y1[0],x1[0]), dst=(y1[1],x1[1]))
    c2 = profile_line(I, src=(y2[0],x2[0]), dst=(y2[1],x2[1]))
    c3 = profile_line(I, src=(y3[0],x3[0]), dst=(y3[1],x3[1]))
    c4 = profile_line(I, src=(y4[0],x4[0]), dst=(y4[1],x4[1]))

    l1 = np.shape(c1)[0]
    l2 = np.shape(c2)[0]
    l3 = np.shape(c3)[0]
    l4 = np.shape(c4)[0]

    # Use zero-padding to get the same size for all arrays
    m = np.max([l1, l2, l3, l4])
    c1 = np.pad(c1, (0,m-l1) ,mode='constant')
    c2 = np.pad(c2, (0,m-l2) ,mode='constant')
    c3 = np.pad(c3, (0,m-l3) ,mode='constant')
    c4 = np.pad(c4, (0,m-l4) ,mode='constant')

    ## Feature extraction

    # linear correlation
    if (typeflag['corr'] or typeflag['local'] or typeflag['texture']):
        roh12 = np.correlate(c1,c2)
        roh13 = np.correlate(c1,c3)
        roh14 = np.correlate(c1,c4)
        roh23 = np.correlate(c2,c3)
        roh24 = np.correlate(c4,c2)
        roh34 = np.correlate(c3,c4)

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
        sd1 = np.std(c1, ddof=1)
        sd2 = np.std(c2, ddof=1)
        sd3 = np.std(c3, ddof=1)
        sd4 = np.std(c4, ddof=1)

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
        sd_dc1 = np.std(dc1, ddof=1)
        sd_dc2 = np.std(dc2, ddof=1)
        sd_dc3 = np.std(dc3, ddof=1)
        sd_dc4 = np.std(dc4, ddof=1)

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

    if not (typeflag['local'] or typeflag['texture']):
        if not typeflag['moments']:
            Out = np.hstack([roh12, roh14, roh24, roh13, roh23, roh34])
        elif not typeflag['corr']:
            Out = np.hstack([m1, m2, m3, m4, m1_5, m2_5, m3_5, m4_5])
        else:
            Out = np.hstack([roh12, roh14, roh24, roh13, roh23, roh34,
                m1, m2, m3, m4, m1_5, m2_5, m3_5, m4_5])
    else:
        Out = np.hstack([roh12, roh14, roh24, roh13, roh23, roh34,
            mean1, mean2, mean3, mean4,
            sd1, sd2, sd3, sd4,
            pr1, pr2, pr3, pr4,
            m1, m2, m3, m4, m1_5, m2_5, m3_5, m4_5,
            p1, p2, p3, p4,
            fft1[:10], fft2[:10], fft3[:10], fft4[:10],
            mean_dc1, mean_dc2, mean_dc3, mean_dc4,
            sd_dc1, sd_dc2, sd_dc3, sd_dc4,
            p_dc2, p_dc1, p_dc3, p_dc4,
            fft_dc1[:10], fft_dc2[:10], fft_dc3[:10], fft_dc4[:10]])

    return Out
