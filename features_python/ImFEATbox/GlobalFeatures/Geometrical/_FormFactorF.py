from skimage import measure
import numpy as np
from __helperCommands import grayscale2bw

def FormFactorF(I, typeflag=None, returnShape=False):
"""
 Input:     - I: A 2D image
            - typeflag: Struct of logicals to permit extracting features
              based on desired characteristics:
                   + typeflag['global']: all features
                   + typeflag['form']: all features
                   + typeflag['corr']: only features based on correlation
              default: all features are being extracted
              For more information see README.txt


 Output:    - Out: A (1x32) vector containing 32 metrics based on the form
                   of detected objects in am image
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


## extract properties of regions in the image

    #if ~exist('typeflag','var')
    if typeflag == None:
        typeflag = dict()
        typeflag['global'] = True
        typeflag['form'] = True
        typeflag['corr'] = True


    if typeflag['global'] or typeflag['form']:
        typeflag['corr'] = True;

    if returnShape:
        if N > 0:
            if (typeflag['global'] or typeflag['form']):
                return (32,1)
            else:
                return (5,1)

        else:
            return (32,1)


    # convert image, image must be real valued!
    BW = grayscale2bw(I)
    L, N = measure.label(BW, background=0, return_num=True)

    # check if any objects have been detected
    if N > 0:
        # extract region properties
        #s = regionprops(L, 'area', 'perimeter', ...
        #    'Orientation', 'Eccentricity', 'EquivDiameter', 'Solidity')

        # cache=True: faster but needs more memory

        # TODO only extract needed properties
        s = measure.regionprops(L, intensity_image=None, cache=True)



        # initialize variables
        roundness = np.zeros(N)
        form = np.zeros(N)
        area_obj = np.zeros(N)
        perimeter_obj = np.zeros(N)
        eccentricity = np.zeros(N)
        orientation = np.zeros(N)
        solidity = np.zeros(N)

        for i in range(N):
            area_obj[i] = s[i].area
            perimeter_obj[i] = s[i].perimeter
            form[i] = 4*np.pi*area_obj[i]/(perimeter_obj[i]^2)
            roundness[i] = (4*s[i].area )/(s[i].equivalent_diameter*np.pi)
            eccentricity[i] = s[i].eccentricity
            orientation[i] = s[i].orientation
            solidity[i] = s[i].solidity

        # avoid Inf values which cause problems in further calculations
        #form(form == Inf) = max(form(form ~= Inf))
        form[np.isinf(form)] = form[~np.isinf(form)].max()

        ## feature extraction

        if (typeflag['global'] or typeflag['form']):
            # mean
            mean_roundness = np.mean(roundness)
            mean_area = np.mean(area_obj)
            mean_perimeter = np.mean(perimeter_obj)
            mean_eccentricity = np.mean(eccentricity)
            mean_orientation = np.mean(orientation)
            mean_solidity = np.mean(solidity)
            mean_form = np.mean(form)

            # standard deviation
            std_roundness = np.std(roundness, ddof=1)
            std_area = np.std(area_obj, ddof=1)
            std_perimeter = np.std(perimeter_obj, ddof=1)
            std_eccentricity = np.std(eccentricity, ddof=1)
            std_orientation = np.std(orientation, ddof=1)
            std_solidity = np.std(solidity, ddof=1)
            std_form = np.std(form, ddof=1)

            # maximum/minimum values (orientation is neglected for this measure)
            max_roundness = np.max(roundness)
            max_area = np.max(area_obj)
            max_perimeter = np.max(perimeter_obj)
            max_eccentricity = np.max(eccentricity)
            max_solidity = np.max(solidity)
            max_form = np.max(form)

            min_roundness = np.min(roundness)
            min_area = np.min(area_obj)
            min_perimeter = np.min(perimeter_obj)
            min_eccentricity = np.min(eccentricity)
            min_solidity = np.min(solidity)
            min_form = np.min(form)

        # calculate features from the correlation between form factor and roundness
        if N > 1:
            #corr_FR = xcorr(form, roundness, 'coeff')

            form_normalized = (form - np.mean(form)) / (np.std(form, ddof=1) * len(form))
            roundness_normalized = (roundness - np.mean(roundness)) / (np.std(roundness) * len(roundness))
            corr_FR = np.correlate(form_normalized,roundness_normalized, "full")

            mean_corr = np.mean(corr_FR)
            std_corr = np.std(corr_FR, ddof=1)
            max_corr = np.max(corr_FR)
            min_corr = np.min(corr_FR)

            # determine cross correlation coefficient
            corrcoef_FR = np.corrcoef(form, roundness)[0,1]
        else:
            mean_corr = 0
            std_corr = 0
            max_corr = 0
            min_corr = 0
            corrcoef_FR = 0



        ## return feature vector
        if (typeflag['global'] or typeflag['form']):
            return np.array([N, mean_eccentricity, mean_solidity, mean_orientation,
                mean_roundness, mean_area, mean_perimeter, mean_form,
                std_eccentricity, std_solidity, std_orientation,
                std_roundness, std_area, std_perimeter, std_form,
                max_eccentricity, max_solidity, max_roundness, max_area, max_perimeter, max_form,
                min_eccentricity, min_solidity, min_roundness, min_area, min_perimeter, min_form,
                corrcoef_FR, mean_corr, std_corr, max_corr, min_corr])
        else:
            return np.array([corrcoef_FR, mean_corr, std_corr, max_corr, min_corr])

    else:
        return np.zeros(32)
