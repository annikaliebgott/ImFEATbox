import numpy as np
from ImFEATbox.__helperCommands import conv2float
from scipy.signal import convolve2d
from scipy.stats import moment

def LawF(I, returnShape=False):
    """
 Input:     - I: A 2D image


 Output:    - Out: A (1x58) vector containing 58 metrics calculated from
                   the 2nd and 4th moments of the results of the image
                   filtered by 25 different filters
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
# Law features are constructed from a set of five one-dimensional filters:
# edge, spot, ripple, wave and lowpass. By applying a one-dimensional
# filter in the horizontal direction followed by a one-dimensional filter
# in the vertical direction, this results in 25 different two-dimensional
# filters.

    if returnShape:
        return (58,1)

    ## filter image horizontal and vertical

    #converte to double
    I = conv2float(I)

    # define 5 filters
    edge = np.array([[-1, -2, 0, 2, 1]])
    low_pass = np.array([[1, 4, 6, 4, 1]])
    spots = np.array([[-1, 0, 2, 0, 1]])
    ripples = np.array([[1, -4, 6, -4, 1]])
    waves = np.array([[-1, 2, 0, -2, 1]])

    # edge filter
    ee = convolve2d_image(edge, edge, I)
    le = convolve2d_image(edge, low_pass, I)
    se = convolve2d_image(edge, spots, I)
    re = convolve2d_image(edge, ripples, I)
    we = convolve2d_image(edge, waves, I)

    #spots
    es = convolve2d_image(spots, edge, I)
    ls = convolve2d_image(spots, low_pass, I)
    ss = convolve2d_image(spots, spots, I)
    rs = convolve2d_image(spots, ripples, I)
    ws = convolve2d_image(spots, waves, I)

    #ripples
    er = convolve2d_image(ripples, edge, I)
    lr = convolve2d_image(ripples, low_pass, I)
    sr = convolve2d_image(ripples, spots, I)
    rr = convolve2d_image(ripples, ripples, I)
    wr = convolve2d_image(ripples, waves, I)

    #waves
    ew = convolve2d_image(waves, edge, I)
    lw = convolve2d_image(waves, low_pass, I)
    sw = convolve2d_image(waves, spots, I)
    rw = convolve2d_image(waves, ripples, I)
    ww = convolve2d_image(waves, waves, I)

    #low pass
    el = convolve2d_image(low_pass, edge, I)
    ll = convolve2d_image(low_pass, low_pass, I)
    sl = convolve2d_image(low_pass, spots, I)
    rl = convolve2d_image(low_pass, ripples, I)
    wl = convolve2d_image(low_pass, waves, I)

    ## feature extraction

    # determine the 2nd and 4th moment of the filtered images
    m2 = np.array([moment(ee,2), moment(se,2), moment(re,2), moment(we,2), moment(le,2),
        moment(es,2), moment(ss,2), moment(rs,2), moment(ws,2), moment(ls,2),
        moment(er,2), moment(sr,2), moment(rr,2), moment(wr,2), moment(lr,2),
        moment(ew,2), moment(sw,2), moment(rw,2), moment(ww,2), moment(lw,2),
        moment(el,2), moment(sl,2), moment(rl,2), moment(wl,2), moment(ll,2)])

    m4 = np.array([moment(ee,4), moment(se,4), moment(re,4), moment(we,4), moment(le,4),
        moment(es,4), moment(ss,4), moment(rs,4), moment(ws,4), moment(ls,4),
        moment(er,4), moment(sr,4), moment(rr,4), moment(wr,4), moment(lr,4),
        moment(ew,4), moment(sw,4), moment(rw,4), moment(ww,4), moment(lw,4),
        moment(el,4), moment(sl,4), moment(rl,4), moment(wl,4), moment(ll,4)])

    # mean of the moments
    mean_m2 = np.mean(m2)
    mean_m4 = np.mean(m4)

    # standard deviation of the moments
    std_m2 = np.std(m2)
    std_m4 = np.std(m4)

    # maximum and minimum values of the moments
    max_m2 = np.max(m2)
    max_m4 = np.max(m4)
    min_m2 = np.min(m2)
    min_m4 = np.min(m4)

    ## return feature vector
    Out = np.array([m2, m4, mean_m2, mean_m4, std_m2, std_m4, max_m2, max_m4, min_m2, min_m4])
    return Out

def convolve2d_image(h1, h2, I):
    # first convolves each column of A with the vector h1
    # and then convolves each row of the result with the vector h2.
    I = convolve2d(I, h1.T, 'same')
    I = convolve2d(I, h2, 'same')
    return I
