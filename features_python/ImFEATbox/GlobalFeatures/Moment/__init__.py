# -*- coding: utf-8 -*-

import ImFEATbox.FeatureGroup as __fg
from ImFEATbox.__Features__ import __Feature

# we use the __ to hide objects from the python IDE,
# you will only see the needed objects

from _ZernikeF import ZernikeF as __ZernikeF
from _HuF import HuF as __HuF
from _AffineMomentsF import AffineMomentsF as __AffineMomentsF

# Moment features
zernike = __Feature('zernike', __ZernikeF, [__fg.global_, __fg.moments])
hu = __Feature('hu', __HuF, [__fg.global_, __fg.moments])
affinemoments = __Feature('affinemoments', __AffineMomentsF, [__fg.global_, __fg.moments])
