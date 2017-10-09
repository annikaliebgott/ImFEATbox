# -*- coding: utf-8 -*-

import ImFEATbox.FeatureGroup as __fg
from ImFEATbox.__Features__ import __Feature

# we use the __ to hide objects from the python IDE,
# you will only see the needed objects

from _FormFactorF import FormFactorF as __FormFactorF


# Geometrical features
formfactor = __Feature('formfactor', __FormFactorF, [__fg.global_, __fg.corr, __fg.form])
