#!/bin/python
# -*- coding: utf-8 -*-

import os
from glob import glob
import ImFEATbox
import itertools
from subprocess import Popen, PIPE, STDOUT


featureWhiteList = []
featureWhiteList.append("svd")


featureFolderList = ["GlobalFeatures", "LocalFeatures"]

testedFeatures = 0
passedFeatures = 0
testedCases = 0

pyFileList = []
for folder in featureFolderList:
    PATH = "ImFEATbox/" + folder
    pl = [os.path.join(r, fn) for r, ds, fs in os.walk(PATH) for fn in fs if fn.endswith('.py')]
    pyFileList += pl

#
# if len(featureWhiteList) > 0:
#     # while testing few features we can use a whitelist
#
#     for keyword in featureWhiteList:
#         for pyFile in pyFileList:
#             pyFileName = pyFile.split(os.sep)[-1]
#             print(keyword)
#             print(pyFileName)
#             if keyword.lower() in pyFileName.lower():
#                 print(keyword.lower() + " in " + pyFile.split(os.sep)[-1].lower())
#                 pass
#             else:
#                 pyFileList.remove(pyFile)
#                 print("remove" + str(pyFile))
#
#
#     print("Using Whitelist: " + str(len(pyFileList)) + " features selected.")
#     print(pyFileList)

report = ""
reportx = ""
report += "######################################################" + os.linesep
report += "############ python ImFEATbox test report ############" + os.linesep
report += "######################################################" + os.linesep + os.linesep
stdout = ""
rmList = set()
for r in pyFileList:
    pyFileName = r.split(os.sep)[-1]
    if pyFileName == "__init__.py" or pyFileName[0] != "_":
        rmList.add(r)

    if len(featureWhiteList) > 0:
        found = False
        for k in featureWhiteList:
            if k.lower() in pyFileName.lower():
                found = True
        if not found:
            rmList.add(r)

#rmList = list(rmList)
for r in rmList:
    pyFileList.remove(r)
print(pyFileList)
print("")
# find corresponding matlab files:
mFileList = []
for r in pyFileList:
    #r = pyFileList[0]
    featurePassed = True
    folder = r[:r.rfind(os.sep)+1]
    pyFileName = r.split(os.sep)[-1]
    imFeature = folder.replace(os.sep, ".") + pyFileName[1:-4].lower()
    reportx += "##### " + str(imFeature) + " #####" + os.linesep

    testedFeatures += 1
    print(imFeature)
    mFileName = folder + pyFileName[1:-2] + "m"
    if not os.path.isfile(mFileName):
        print("Error! File " + mFileName + " does not exist!")

    getPyParameters = "pyParameters = " + imFeature + ".getInputParameters()"
    pyMethod =  imFeature + ".cFeatures"


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
                        l = l.replace(os.linesep, '')
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
    plotParameter = False
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
        elif par == "plotflag":
            plotParameter = True
        else:
            print("unknown parameter: " + par)
            quit()
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

    reportx += "relevant typeFlags:" +os.linesep
    reportx += str(typeFlagList) + os.linesep

    typeFlagBools = list(itertools.product([False, True], repeat=len(typeFlagList)))
    #print(typeFlagBools)





    # now the test cases:
    mMethod = mFileName[:-2].split(os.sep)[-1]




    # single script testing:

    for j in range(0, len(typeFlagBools)):
        testedCases += 1
        parameterSetup = ""
        cmdPyShape = "pyShape = " + imFeature + ".getFeatureShape"
        pythonScript = "import ImFEATbox, csv" + os.linesep
        pythonScript += "import numpy as np" + os.linesep
        pythonScript += "with open('testimg.csv', 'r') as csvfile:" + os.linesep
        pythonScript += "\tI = np.array(list(csv.reader(csvfile, delimiter=','))).astype(np.float)" + os.linesep

        matlabScript = "I = csvread('testimg.csv');" + os.linesep
        if plotParameter:
            pythonScript += "plotflag = False" + os.linesep
            matlabScript += "plotflag = false;" + os.linesep
            parameterSetup += "plotflag = False" + os.linesep
            exec("plotflag = False")
        if gradTypeParameter:
            pythonScript += "gradtype = dict()" + os.linesep
            exec("gradtype = dict()")
        if len(typeFlagList) > 0:
            pythonScript += "typeflag = dict()" + os.linesep
            exec("typeflag = dict()")
        for i in range(0, len(typeFlagList)):
            matlabScript += typeFlagList[i] + " = " + str(typeFlagBools[j][i]).lower() + ";" + os.linesep
            pythonScript += typeFlagList[i].replace(".", "['") + "'] = " + str(typeFlagBools[j][i]) + "" + os.linesep
            parameterSetup += typeFlagList[i].replace(".", "['") + "'] = " + str(typeFlagBools[j][i]) + "" + os.linesep
            exec(typeFlagList[i].replace(".", "['") + "'] = " + str(typeFlagBools[j][i]))
        matlabScript += "addpath('" + folder + "');" + os.linesep
        matlabScript += "Out = " + mMethod + "(I, "
        for par in mParameterList:
            matlabScript += par + ", "
        matlabScript = matlabScript[:-2] + ");" + os.linesep
        pythonScript += "out = " + pyMethod + "(I, "
        cmdPyShape += "(None, "
        for par in pyParameters:
            pythonScript += par + ", "
            cmdPyShape += par + ", "
        pythonScript = pythonScript[:-2] + ")" + os.linesep
        cmdPyShape = cmdPyShape[:-2] + ")"
        matlabScript += "csvwrite('matlab-out.csv', Out);"
        pythonScript += "with open('python-out.csv', 'wb') as csvfile:" + os.linesep
        pythonScript += "\twriter = csv.writer(csvfile, delimiter=',', escapechar=' ', quoting=csv.QUOTE_NONE)" + os.linesep
        pythonScript += "\twriter.writerow(out)"
        with open ("tmp.m", 'w+') as f:
            f.write(matlabScript)
        with open ("tmp.py", 'w+') as f:
            f.write(pythonScript)

        # read shape, stored in variable: pyShape
        #print(cmdPyShape)
        exec(cmdPyShape)

        reportx += "test case #" + str(j+1) + "/" + str(len(typeFlagBools)) + ":" + os.linesep
        reportx += parameterSetup

        # now execute both scripts:
        execfile("tmp.py")
        matlabParameter = " -nojvm -nodisplay -nosplash -r \"try, run('tmp.m'), , catch, exit, end, exit\""
        #os.system("matlab" + matlabParameter)
        mclo = Popen(["matlab", matlabParameter], stdout=PIPE, stderr=STDOUT)
        stdout += mclo.communicate()[0]



        # read both Output csv
        mComplexValues = False
        with open('matlab-out.csv', 'r') as csvfile:
            mRead = list(csv.reader(csvfile, delimiter=','))[0]
            #print("len=" + str(len(mRead)))
            for i in range(0, len(mRead)):
                #print(mRead[i])
                if "i" in mRead[i]:
                    # in python is imaginary symbol j
                    mRead[i] = mRead[i].replace("i","j")
                    mComplexValues = True

            #print("mRead:")
            #print(mRead)
            if mComplexValues:
                mOut = np.array(mRead).astype(np.complex128)
            else:
                mOut = np.array(mRead).astype(np.float)


        pComplexValues = False
        with open('python-out.csv', 'r') as csvfile:
            pRead = list(csv.reader(csvfile, delimiter=','))[0]
            for i in range(0, len(pRead)):
                if "j" in pRead[i]:
                    pComplexValues = True
            if pComplexValues:
                pOut = np.array(pRead).astype(np.complex128)
            else:
                pOut = np.array(pRead).astype(np.float)



        if len(pyShape) == 1:
            pyShapeNum = pyShape[0]
        else:
            # if we use matrix output
            pyShapeNum = pyShape[0]*pyShape[1]

        # check if returnShape is correct:

        if pyShapeNum == len(pOut) == len(mOut):
            reportx += "\t* shape test passed: " + str(pyShapeNum) + os.linesep
        else:
            featurePassed = False
            reportx += "\t* shape test FAILED !!!" + os.linesep
            print("mOut: " + str(len(mOut)) + ", pOut: " + str(len(pOut)) + ", shape: " + str(pyShape))
            print("length wrong!")
            print(parameterSetup)

        diff = pOut - mOut
        diff_norm = np.abs(diff/mOut)
        max_diff_norm = np.max(np.abs(diff_norm))*100

        if max_diff_norm < 1:
            reportx += "\t* value test passed: all value deviations < " + str(max_diff_norm) + " %" + os.linesep
        else:
            featurePassed = False
            reportx += "\t* value test FAILED !!!" + os.linesep
            reportx += "maximum deviation: " + str(max_diff_norm) + " %" + os.linesep
            reportx += os.linesep
            reportx += "Matlab:" + os.linesep
            reportx += str(mOut) + os.linesep
            reportx += "Python:" + os.linesep
            reportx += str(pOut) + os.linesep
            reportx += "(Python-Matlab):" + os.linesep
            reportx += str(diff) + os.linesep
            reportx += "(Python-Matlab)/Matlab:" + os.linesep
            reportx += str(diff_norm) + os.linesep
            diffindex1 = np.where(np.abs(diff_norm) > 0.01)[0]
            diffindex5 = np.where(np.abs(diff_norm) > 0.05)[0]
            diffindex10 = np.where(np.abs(diff_norm) > 0.1)[0]
            reportx += "feature indices with > 1% Error:" + os.linesep
            reportx += str(diffindex1) + os.linesep
            reportx += "feature indices with > 5% Error:" + os.linesep
            reportx += str(diffindex5) + os.linesep
            reportx += "feature indices with > 10% Error:" + os.linesep
            reportx += str(diffindex10) + os.linesep
        reportx += os.linesep
    if featurePassed:
        passedFeatures += 1

report += "Overall Stats:" + os.linesep
#report += "\t* " + str(testedFeatures) + " features tested." + os.linesep
report += "\t* " + str(testedCases) + " cases tested." + os.linesep
report += "\t* " + str(passedFeatures) + "/" + str(testedFeatures) + " features passed." + os.linesep
report += os.linesep
report += reportx

with open('matlab.log', 'w+') as matlablog:
    matlablog.write(stdout)
with open('testreport.txt', 'w+') as testreportFile:
    testreportFile.write(report)

print(os.linesep + os.linesep + "***** done *****" + os.linesep)
