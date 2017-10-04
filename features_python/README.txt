Work in progress

Source Matlab files are from 9.04.2017
Matlab files still remain in the folder until
updated to current version and properly tested


************************************************************************
This file is part of ImFEATbox, a toolbox for image feature extraction and
analysis. Implemented by the Department of Diagnostic and Interventional
Radiology, University Hospital of Tuebingen, Germany and the Institute of
Signal Processing and System Theory University of Stuttgart, Germany.
Last modified: February 2017

Available online at:
https://github.com/annikaliebgott/ImFEATbox

Contact: annika.liebgott@iss.uni-stuttgart.de
		 thomas.kuestner@iss.uni-stuttgart.de
************************************************************************

Quick Tutorial:

To be able to import this library, either have it as a subdirectory or add it
to python path.
In future there will be a PIP package which does all that automatically.

>>> import ImFEATbox

Now yow can access the implemented features and functions in a hierarchic structure:

>>> ImFEATbox.
ImFEATbox.FeatureGroup     ImFEATbox.getFeatureNames
ImFEATbox.GlobalFeatures   ImFEATbox.getFeatures
ImFEATbox.LocalFeatures

To see the doc string and get help for methods use the '?' after a method:
>>> ImFEATbox.getFeatures?
Type:        function
String form: <function getFeatures at 0x7f927b9d3050>
File:        /home/heiko/git/ImFEATbox/features_python/ImFEATbox/__Features__.py
Definition:  ImFEATbox.getFeatures(nameList=None)
Docstring:
Input:  - nameList: A list of feature names
            e.g. items of getFeatureNames()
          Default: None
            -> returns all available/implemented features

Output:    - A list of __Feature objects

Now we load a test image and extract all available features:
(image loader function will be implemented in this library in the future)

>>> import csv
>>> import numpy as np
>>> with open('testimg.csv', 'r') as csvfile:
>>> 	I = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)

Iterate through all available features and store output in a list:

>>> Out = []
>>> featureList = ImFEATbox.getFeatures()
>>> for feature in featureList:
>>>     Out.append(feature.cFeatures(I))

Now flatten the list and convert it to numpy array:

>>> Out = np.hstack(Out)

We have now a big feature vector for that one image with length of 1170:
>>> Out.shape
(1170,)

For machine learning we want a dataset of X (training vectors) and Y (class labels)
Lets assume we have a set of images labeled as cancer and one for non-cancer:

>>> cancerList = loadImages("cancer/*") # pseudo-code
>>> noncancerList = loadImages("noncancer/*") # pseudo-code
>>> imgList = noncancerList + cancerList # appending
>>> Y = [0]*len(noncancerList) + [1]*len(cancerList)
>>> X = []
>>> for img in imgList:
>>> 	Out = []
>>> 	for feature in featureList:
>>> 		Out.append(feature.cFeatures(I))
>>> 	X.append(np.hstack(Out))

Now we have our X feature set and corresponding Y labels
X[0] is labeled with Y[0] and so on... X[i] belongs to class Y[i]

Be aware that the feature vectors are not normalized yet.
