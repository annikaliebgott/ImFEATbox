# -*- coding: utf-8 -*-
import numpy as np

__float_dtype = np.float64
__complex_dtype = np.complex128


def _complex_dtype():
    return __complex_dtype

def _float_dtype():
    return __float_dtype

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
    return np.array(Image, dtype='float64')

def isColorImage(Image):
    """
    returns True if the input parameter is an Image in color.
    otherwise returns False.
    working for 3D and 2D images

    checking by the shape
    """
    if len(np.shape(Image)) == 4 and np.shape(Image)[3] == 3: #image is 3D and color
        return True
    elif len(np.shape(Image)) == 3 and np.shape(Image)[2] == 3: #image is 2D and color
        return True
    else:
        return False
