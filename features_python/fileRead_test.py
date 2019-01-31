#!/bin/python
# -*- coding: utf-8 -*-

import os
#from glob import glob
import ImFEATbox
import numpy as np
#np.set_printoptions(threshold=np.nan,precision=2)

testImageFolder = "/home/heiko/Nextcloud/data/python/imfeatbox_data"
testImage = "/home/heiko/Nextcloud/data/python/imfeatbox_data/testimage.IMA"

il = ImFEATbox.ImageLoader()

il.loadFolder(testImageFolder)
