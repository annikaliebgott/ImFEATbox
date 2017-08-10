# -*- coding: utf-8 -*-

import numpy as np

def rgb2grayscale(Image):
    """
    converts an rgb image to grayscale using the matlab weighting:
    0.2989 * R + 0.5870 * G + 0.1140 * B
    """
    if type(Image) != np.ndarray:
        Image = np.array(Image)
    return np.dot(Image[...,:3], [0.2989, 0.5870, 0.1140])

def conv2float(Image):
    """
    converts any array to a numpy array with float64 datatype
    """
    return np.array(I, dtype='float64')
