# -*- coding: utf-8 -*-

import Line
import Point

import ImFEATbox.FeatureGroup as __fg
from ImFEATbox.__Features__ import __Feature

# we use the __ to hide objects from the python IDE,
# you will only see the needed objects

from _HarrisF import HarrisF as __HarrisF

# Point features
harris = __Feature('harris', __HarrisF, [__fg.local, __fg.texture])
