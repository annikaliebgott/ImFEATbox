# -*- coding: utf-8 -*-

import ImFEATbox.FeatureGroup as __fg
from ImFEATbox.__Features__ import __Feature

# we use the __ to hide objects from the python IDE,
# you will only see the needed objects

from _IntensityF import IntensityF as __IntensityF
from _GradientF import GradientF as __GradientF
from _SVDF import SVDF as __SVDF
from _HistogramF import HistogramF as __HistogramF

# Intensity features
intensity = __Feature('intensity', __IntensityF, [__fg.global_, __fg.corr, __fg.texture, __fg.entropy])
histogram = __Feature('gradient', __HistogramF, [__fg.global_, __fg.texture, __fg.entropy])
svd = __Feature('svd', __SVDF, [__fg.global_, __fg.texture])
gradient = __Feature('gradient', __GradientF, [__fg.global_, __fg.gradient, __fg.texture, __fg.entropy])
