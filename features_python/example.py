#!/bin/python
# -*- coding: utf-8 -*-

# Packages:
# sudo apt-get install python-scipy
# sudo apt-get install python-skimage



#import ImFEATbox
#from PIL import Image
import Image
import numpy as np
import ImFEATbox
from ImFEATbox.__helperCommands import rgb2grayscale

import matplotlib.pyplot as plt

#print(ImFEATbox.getFeatureNames())


#img = Image.open("n.jpg")
I = Image.open("lena.tiff")

I = rgb2grayscale(I)

f=ImFEATbox.GlobalFeatures.Intensity.gradient.cFeatures()




shape = np.shape(I)
print(shape)




plt.imshow(I, cmap='gray')
plt.show()
