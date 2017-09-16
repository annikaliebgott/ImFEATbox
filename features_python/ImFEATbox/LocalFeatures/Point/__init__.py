# -*- coding: utf-8 -*-

import ImFEATbox.FeatureGroup as __fg
from ImFEATbox.__Features__ import __Feature

# we use the __ to hide objects from the python IDE,
# you will only see the needed objects

from _LawF import LawF as __LawF

# Point features
law = __Feature('LawF', __LawF, [__fg.local, __fg.moments, __fg.texture])
