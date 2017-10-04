import numpy as np

def RunLengthF(I):
"""
     Input:     - I: A 2D image


     Output:    - Out: A (1x44) vector containing 44 metrics calculated based
                  on run length of image structures

     Note: This function calls grayrlmatrix.py and grayrlprops.py, which are
           also included in ImFEATbox and have been implemented by:
           Xunkai Wei <xunkai.wei@gmail.com>, Beijing Aeronautical Technology
           Research Center
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

    # grayrlmatrix.m can not process complex values
    I = double(np.real(I))

    [GLRLMS,~] = grayrlmatrix(I,NumLevels=255, GrayLimits=[np.min(I[:]), np.max(I[:])])
    VectorDegree = grayrlprops(GLRLMS)
    Out = [VectorDegree[1,:], VectorDegree[2,:], VectorDegree[3,:], VectorDegree[4,:]]

    return Out
