#!/bin/python
# -*- coding: utf-8 -*-

import os
from glob import glob
import ImFEATbox
import itertools

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
    #r = pyFileList[0]

    folder = r[:r.rfind(os.sep)+1]
    pyFileName = r.split(os.sep)[-1]
    imFeature = folder.replace(os.sep, ".") + pyFileName[1:-4].lower()
    print(imFeature)
    mFileName = folder + pyFileName[1:-2] + "m"
    if not os.path.isfile(mFileName):
        print("Error! File " + mFileName + " does not exist!")

    getPyParameters = "pyParameters = " + imFeature + ".getInputParameters()"
    pyMethod =  imFeature + ".cFeatures"
    cmdPyShape = "pyShape = " + imFeature + ".getFeatureShape()"
    exec(cmdPyShape)
    exec(getPyParameters)
    #print(pyParameters)
    #for par in pyParameters:

    mFullFileName = os.path.abspath(mFileName)
    #print(mFullFileName)
    mFileList.append(mFileName)
    #print(os.path.abspath(os.path.curdir))

    # grab matlab parameters
    typeFlagList = set()
    mParameterList = []
    with open (mFileName, 'r') as f:
        mHeadfound = False
        for line in f:
            lineSplit = line.split(" ")
            if "typeflag." in line and "if" in line:
                for l in lineSplit:
                    if "typeflag." in l:
                        l = l.replace('(', '')
                        l = l.replace(')', '')
                        l = l.replace('~', '')
                        l = l.replace('\n', '')
                        typeFlagList.add(l)
            if not mHeadfound:
                if "=" in line and "function" in line and "(" in line and ")" in line:
                    mHeadfound = True
                    i1 = line.find("(")
                    i2 = line.find(")")
                    pars = line[i1+1:i2].replace(" ", "")
                    mParameterList = pars.split(",")
                    # should be the method head



    # looking for unknown parameters:
    gradTypeParameter = False
    for par in mParameterList:
        if par.lower() == "i" or par.lower() == "image":
            pass
        elif par == "typeflag":
            pass
        elif par == "gradtype":
            # ok keep it and since it's also boolean, add it to typeFlagList
            # gradtype.first, gradtype.second
            typeFlagList.add("gradtype.first")
            typeFlagList.add("gradtype.second")
            gradTypeParameter = True
        else:
            print("unknown parameter: " + par)
    #print(typeFlagList)

    # now convert to python-parameters:

    if "returnShape" in pyParameters:
        pyParameters.remove("returnShape")
    else:
        print("error! forget returnShape? " + imFeature)
        quit()
    if "I" in pyParameters:
        pyParameters.remove("I")
    if "I" in mParameterList:
        mParameterList.remove("I")

    if len(mParameterList) > len(pyParameters):
        print("error: missing paramters! " + imFeature)
        print("python: " + str(pyParameters))
        print("matlab: " + str(mParameterList))
        quit()


    # create test cases:
    typeFlagList = list(typeFlagList)

    typeFlagBools = list(itertools.combinations_with_replacement([True, False], len(typeFlagList)))






    # now the test cases:
    mMethod = mFileName[:-2]




    # single script testing:

    for j in range(0, len(typeFlagBools)):
        pythonScript = "import ImFEATbox, csv\n"
        pythonScript += "import numpy as np\n"
        pythonScript += "with open('testimg.csv', 'r') as csvfile:\n"
        pythonScript += "\tI = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)\n"

        matlabScript = "I = csvread('testimg.csv');\n"
        if gradTypeParameter:
            pythonScript += "gradtype = dict()\n"
        if len(typeFlagList) > 0:
            pythonScript += "typeflag = dict()\n"
        for i in range(0, len(typeFlagList)):
            matlabScript += typeFlagList[i] + " = " + str(typeFlagBools[j][i]).lower() + ";\n"
            pythonScript += typeFlagList[i].replace(".", "['") + "'] = " + str(typeFlagBools[j][i]) + "\n"
        matlabScript += "addpath('" + folder + "');\n"
        matlabScript += "Out = " + mMethod + "(I, "
        for par in mParameterList:
            matlabScript += par + ", "
        matlabScript = matlabScript[:-2] + ");\n"
        pythonScript += "out = " + pyMethod + "(I, "
        for par in pyParameters:
            pythonScript += par + ", "
        pythonScript = pythonScript[:-2] + ");\n"
        matlabScript += "csvwrite('matlab-out.csv', Out);"
        pythonScript += "with open('python-out.csv', 'wb') as csvfile:\n"
        pythonScript += "\twriter = csv.writer(csvfile, delimiter=',')\n"
        pythonScript += "\twriter.writerow(out)"
        with open ("tmp.m", 'w+') as f:
            f.write(matlabScript)
        with open ("tmp.py", 'w+') as f:
            f.write(pythonScript)

        # now execute both scripts:
        execfile("tmp.py")

        # check if returnShape is correct:

        if len(pyShape) == 1 == 1:
            pass
