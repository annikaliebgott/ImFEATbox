#!/bin/python
# -*- coding: utf-8 -*-

#import ImFEATbox
#from PIL import Image
import Image
import numpy as np
import ImFEATbox
from ImFEATbox.__helperCommands import rgb2grayscale
import csv
import matplotlib.pyplot as plt

#print(ImFEATbox.getFeatureNames())

# load test image
with open('testimg.csv', 'r') as csvfile:
    I = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)

#I = rgb2grayscale(I)

out_python = ImFEATbox.GlobalFeatures.Intensity.gradient.cFeatures(I)


###############################################################################
# * Matlab Code for CSV feature extract:
# csvwrite('matlab-out.csv', Out)
###############################################################################

# .. now read matlab csv:
with open('matlab-out.csv', 'r') as csvfile:
    out_matlab = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float).ravel()

# now compare matlab and python output:

if len(out_python) != len(out_matlab):
    print("Problem: not same # of features: python: " + str(len(out_python)) + ", matlab: " + str(len(out_matlab)))
    print("quit")
    quit()


print("shape matlab: " + str(np.shape(out_matlab)) + ", shape python: " + str(np.shape(out_python)))

# we see matlab as reference code
diff = out_matlab - out_python
maxval_matlab = np.max(out_matlab)
maxval_python = np.max(out_python)
minval_matlab = np.min(out_matlab)
minval_python = np.min(out_python)
valrange_matlab = maxval_matlab - minval_matlab
valrange_python = maxval_python - minval_python

#print(np.abs(diff))




# max value differs 5% of value range
if abs(maxval_matlab - maxval_python) > 0.05 * valrange_matlab:
    print("Problem: maximum differs > 5 percent of value range")

# corresponding indices
diffindex01 = np.where(np.abs(diff) < 0.001*valrange_matlab)[0]
diffindex1 = np.where(np.abs(diff) > 0.01*valrange_matlab)[0]
diffindex5 = np.where(np.abs(diff) > 0.05*valrange_matlab)[0]
diffindex10 = np.where(np.abs(diff) > 0.1*valrange_matlab)[0]
diffindex = np.where(np.abs(diff) > 0)[0]

diffpercentage = (diff/valrange_matlab)*100

maxdiffpercentage = np.max(np.abs(diffpercentage))

print("Result:")
print("\t* " + str(len(diffindex01)) + " values match < 0.1 percent")
print("\t* maximum difference: " + str(maxdiffpercentage) + " percent")
print("\t* " + str(len(diffindex1)) + " values differ > 1 percent")
print("\t* " + str(len(diffindex5)) + " values differ > 5 percent")
print("\t* " + str(len(diffindex10)) + " values differ > 10 percent")

#print(diffindex)
#for i in diffindex[0]:
#    print("i=" + str(i) + " : " + str(dif[i]))


#plt.imshow(I, cmap='gray')

plt.bar(range(1,len(diff)+1), diffpercentage)
plt.show()
