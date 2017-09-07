#!/bin/python
# -*- coding: utf-8 -*-

import os
from glob import glob

featureFolderList = ["GlobalFeatures", "LocalFeatures"]
PATH = "ImFEATbox/GlobalFeatures"
#result = [y for x in os.walk(PATH) for y in glob(os.path.join(x[0], '*.py'))]

pyFileList = [os.path.join(r, fn) for r, ds, fs in os.walk(PATH) for fn in fs if fn.endswith('.py')]


rmList = []
for r in pyFileList:
    pyFileName = r.split(os.sep)[-1]
    if pyFileName == "__init__.py" or pyFileName[0] != "_":
        rmList.append(r)
for r in rmList:
    pyFileList.remove(r)

# find corresponding matlab files:
mFileList = []
for r in pyFileList:
    pyFileName = r.split(os.sep)[-1]
    folder = r[:r.rfind(os.sep)+1]
    mFileName = folder + pyFileName[1:-2] + "m"
    if not os.path.isfile(mFileName):
        print("Error! File " + mFileName + " does not exist!")
    mFullFileName = os.path.abspath(mFileName)
    print(mFullFileName)
    mFileList.append(mFileName)
    print(os.path.abspath(os.path.curdir))
