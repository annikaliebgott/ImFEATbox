# -*- coding: utf-8 -*-

import ImFEATbox.FeatureGroup as __fg
from ImFEATbox.__Features__ import __Feature

# we use the __ to hide objects from the python IDE,
# you will only see the needed objects

from _LineProfileF import LineProfileF as __LineProfileF

# Line features
lineprofile = __Feature('lineProfile', __LineProfileF, [__fg.local, __fg.corr, __fg.moments, __fg.texture])
