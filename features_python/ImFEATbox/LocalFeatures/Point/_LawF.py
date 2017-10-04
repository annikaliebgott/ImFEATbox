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
    m2 = np.hstack([moment(ee.ravel(),2), moment(se.ravel(),2), moment(re.ravel(),2), moment(we.ravel(),2), moment(le.ravel(),2),
        moment(es.ravel(),2), moment(ss.ravel(),2), moment(rs.ravel(),2), moment(ws.ravel(),2), moment(ls.ravel(),2),
        moment(er.ravel(),2), moment(sr.ravel(),2), moment(rr.ravel(),2), moment(wr.ravel(),2), moment(lr.ravel(),2),
        moment(ew.ravel(),2), moment(sw.ravel(),2), moment(rw.ravel(),2), moment(ww.ravel(),2), moment(lw.ravel(),2),
        moment(el.ravel(),2), moment(sl.ravel(),2), moment(rl.ravel(),2), moment(wl.ravel(),2), moment(ll.ravel(),2)])

    m4 = np.hstack([moment(ee.ravel(),4), moment(se.ravel(),4), moment(re.ravel(),4), moment(we.ravel(),4), moment(le.ravel(),4),
        moment(es.ravel(),4), moment(ss.ravel(),4), moment(rs.ravel(),4), moment(ws.ravel(),4), moment(ls.ravel(),4),
        moment(er.ravel(),4), moment(sr.ravel(),4), moment(rr.ravel(),4), moment(wr.ravel(),4), moment(lr.ravel(),4),
        moment(ew.ravel(),4), moment(sw.ravel(),4), moment(rw.ravel(),4), moment(ww.ravel(),4), moment(lw.ravel(),4),
        moment(el.ravel(),4), moment(sl.ravel(),4), moment(rl.ravel(),4), moment(wl.ravel(),4), moment(ll.ravel(),4)])

    # mean of the moments
    mean_m2 = np.mean(m2)
    mean_m4 = np.mean(m4)

    # standard deviation of the moments
    std_m2 = np.std(m2, ddof=1)
    std_m4 = np.std(m4, ddof=1)

    # maximum and minimum values of the moments
    max_m2 = np.max(m2)
    max_m4 = np.max(m4)
    min_m2 = np.min(m2)
    min_m4 = np.min(m4)

    ## return feature vector
    Out = np.hstack([m2, m4, mean_m2, mean_m4, std_m2, std_m4, max_m2, max_m4, min_m2, min_m4])
    return Out

def convolve2d_image(h1, h2, I):
    # first convolves each column of A with the vector h1
    # and then convolves each row of the result with the vector h2.
    I = convolve2d(I, h1.T, boundary='fill', mode='same')
    I = convolve2d(I, h2, boundary='fill', mode='same')
    return I
