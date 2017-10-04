# -*- coding: utf-8 -*-

#from ImFEATbox.FeatureGroup import _global as global

#from ImFEATbox.__Features__ import __grouptypes, __FeatureList, getFeatures, getFeatureNames

# this file is to register features to feature groups.
# a group can have multiple features -> featureList
# a feature can also have multiple groups


__groupList = []

class __Group(object):
    """
    name: [String] name of the feature
    """

    def __init__(self, groupList):
        self.featureList = set() # unique items
        groupList.append(self)

    def addFeature(self, feature):
        self.featureList.add(feature)

    def getFeatures(self):
        return self.featureList

# create grouptypes
# a list of all grouptypes shall be in __groupList

undefined = __Group(__groupList)
global_ = __Group(__groupList)
local = __Group(__groupList)
corr = __Group(__groupList)
gradient = __Group(__groupList)
moments = __Group(__groupList)
texture = __Group(__groupList)
form = __Group(__groupList)
entropy = __Group(__groupList)
transform = __Group(__groupList)
