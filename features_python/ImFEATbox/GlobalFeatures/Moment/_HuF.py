import numpy as np
import math
from ImFEATbox.__helperCommands import conv2float

def HuF(I):
    """
     Input:     - I: A 2D image


     Output:    - Out: A (1x8) vector containing the first 8 Hu moments of
                       order 3
    """
    # ************************************************************************
    # Modified for MRI feature extraction by the Department of Diagnostic
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
    # Implementation by:    Ishrat Badami, Computer Graphics Department,
    #                       University of Bonn

    I = np.array(I, dtype='float')
    height, width = np.shape(I)

    # define a coordinate system for the image
    xgrid = np.repeat(np.array(range(-np.floor(height/2),np.ceil(height/2))-1).T, [1, width])
    ygrid = np.repeat(np.array(range(-np.floor(width/2),np.ceil(width/2))-1),[height,1])

    x_bar, y_bar = centerOfMass(I,xgrid,ygrid)

    # normalize coordinate system by subtracting mean
    xnorm = x_bar - xgrid
    ynorm = y_bar - ygrid

    # central moments
    mu_11 = central_moments(I, xnorm, ynorm, 1, 1)
    mu_20 = central_moments(I, xnorm, ynorm, 2, 0)
    mu_02 = central_moments(I, xnorm, ynorm, 0, 2)
    mu_21 = central_moments(I, xnorm, ynorm, 2, 1)
    mu_12 = central_moments(I, xnorm, ynorm, 1, 2)
    mu_03 = central_moments(I, xnorm, ynorm, 0, 3)
    mu_30 = central_moments(I, xnorm, ynorm, 3, 0)


    #calculate first 8 hu moments of order 3
    I_1 = mu_20 + mu_02
    I_2 = (mu_20 - mu_02)^2 + 4*mu_11
    I_3 = (mu_30 - 3*mu_12)^2 + (mu_03 - 3*mu_21)^2
    I_4 = (mu_30 + mu_12)^2 + (mu_03 + mu_21)^2
    I_5 = (mu_30 - 3*mu_12)*(mu_30 + mu_12)*((mu_30 + mu_12)^2 -3*(mu_21 + mu_03)^2) + (3*mu_21 - mu_03)*(mu_21 + mu_03)*(3*(mu_30 + mu_12)^2 - (mu_03 + mu_21)^2)
    I_6 = (mu_20 - mu_02)*((mu_30 + mu_12)^2 - (mu_21 + mu_03)^2) + 4*(mu_30 + mu_12)*(mu_21 + mu_03)
    I_7 = (3*mu_21 - mu_03)*(mu_30 + mu_12)*((mu_30 + mu_12)^2 - 3*(mu_21 + mu_03)^2) + (mu_30 - 3*mu_12)*(mu_21 + mu_03)*(3*(mu_30 + mu_12)^2 - (mu_03 + mu_21)^2)
    I_8 = mu_11*(mu_30 + mu_12)^2 - (mu_03 + mu_21)^2 - (mu_20 - mu_02)*(mu_30 + mu_12)*(mu_21 + mu_03)

    Out = np.concatenate(I_1, I_2, I_3, I_4, I_5, I_6, I_7, I_8)

    return Out

    # calculate scale invariant central moments
def central_moments(I,xnorm,ynorm,p,q):

    image = conv2float(I)
    cm = np.sum(np.sum(np.power(xnorm, p)*np.power(ynorm, q)*image))
    cm_00 = sum(image[:]) #this is same as mu(0,0)
    # normalise moments for scale invariance
    cm = cm / np.power(cm_00, (1+(p+q)/2))

    return cm

    # calculate center of mass
def centerOfMass(I,xgrid,ygrid):
    # very small constant to prevent dividing by zero
    eps = math.pow(10,-6)

    x_bar = np.sum(np.sum((xgrid*I)))/(np.sum(I[:]+eps))
    y_bar = np.sum(np.sum((ygrid*I)))/(np.sum(I[:]+eps))

    return x_bar, y_bar
