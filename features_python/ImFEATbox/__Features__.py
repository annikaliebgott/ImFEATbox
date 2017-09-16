# -*- coding: utf-8 -*-

import inspect

class __Feature(object):
    """
    name: [String] name of the feature
    cFeatures: method to calculate the features
    """
    _FeatureList = []

    def _getFeatureList(self):
        return type(self)._FeatureList

    def __init__(self, name, handleMethod, group):
        self.name = name
        self.cFeatures = handleMethod
        self.group = group
        self._getFeatureList().append(self)

    def getInputParameters(self):
        return inspect.getargspec(self.cFeatures)[0]

    def getFeatureShape(self, *args):
        """
        returns a tuple (x,y) with the shape of the feature vector.
        """
        #if "typeflag" in self.getInputParameters():
        return self.cFeatures(*args, returnShape=True)
        #else:
        #    return self.cFeatures(None, returnShape=True)

    def getFeatureSize(self, typeflag=None):
        """
        returns the number of elements in the feature vector.
        """

    def __repr__(self):
        """ what the IDE prints """
        return "Feature object: " + self.name


def __findFeaturefromName(featureName):
    for feature in __Feature._FeatureList:
        if featureName == feature.name:
            return feature
    raise ValueError('Unknown feature name: ' + featureName)
    return None

def getFeatures(nameList=None):
    """
    Input:  - nameList: A list of feature names
                e.g. items of getFeatureNames()
              Default: None
                -> returns all available/implemented features

    Output:    - A list of __Feature objects
    """
    if nameList == None:
        return __Feature._FeatureList
    else:
        featureList = []
        for featureName in nameList:
            featureList.append(__findFeaturefromName(featureName))
        return featureList




def getFeatureNames():
    """returns a String list of all available/implemented features"""
    nameList = []
    for f in __Feature._FeatureList:
        nameList.append(f.name)
    return nameList
