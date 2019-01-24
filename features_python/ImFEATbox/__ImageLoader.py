# -*- coding: utf-8 -*-

import numpy as np
import os
import dicom # pip install pydicom

class __ImageFile(object):
    """
    loads image files to a list and add labels
    split
    """



    def __init__(self, fileLocation, label):
        self.fileLocation = fileLocation
        self.label = label

    def getImage(self):
        #TODO
        pass

class __ImageLoader(object):
    """
    loads image files to a list and add labels
    split
    """



    def __init__(self, forceDatatype=None, forceNormalize=None):
        """
        create an Image loader object which treats loaded images
        in the same way. All loaded images e.g. can be normalized to
        8 bit (256) integer values.

        """

        self.forceDatatype = forceDatatype
        self.forceNormalize = forceNormalize
        self.imageList = []



    def normalizeImages(self, dataType):
        """
        ok
        """

        for imgf in imageList:

        # TODO rework this, a complete loaded image set
        # should be processed
        if self.forceNormalize != None:
            numLevels = None
            if self.forceDatatype == np.int8:
                numLevels = 256

            if numLevels == None:
                # scale to self.forceNormalize (may be 1.0)

                pass
            else:
                GrayLimits = [np.min(Image), np.max(Image)]
                slope = (numLevels-1) / float(GrayLimits[1] - GrayLimits[0])
                intercept = 0 - (slope*(GrayLimits[0]))
                SI = np.round((slope*Image+intercept), decimals=0).astype(self.forceDatatype)




    def loadFolder(folder, stringSelect=None, label=None, subFolderLabel=False):
        """
        This method loads image sets.
            - folder: relative or absolute path to image folder
            - stringSelect: may specifying content of the filenames to load
                String: e.g. "*.jpg", "IMG100*"
            - label:    String: specifying a label for all images to load
            - subFolderLabel:
                True:   loading every subfolder and adding labels according to
                        the subfolder name.
        """
        folderConent = os.listdir(folder)
            if subFolderLabel:
                subFolderList = []
                for f in folderConent:
                    if os.path.isdir(f):
                        subFolderList.append(f)
                for sf in subFolderList:
                    label = sf
                    sfLocation = folder + os.sep + sf
                    for sbf in os.listdir(sfLocation):
                        if sbf meets Stringselect and is file: # TODO
                            iLocation = sfLocation + os.sep + sbf
                            self.imageList.append(__ImageFile(iLocation, label))

            elif labels != None: # no subFolders
                for f in folderConent:
                    if f meets Stringselect and is file: # TODO
                        iLocation = folder + os.sep + f
                        self.imageList.append(__ImageFile(iLocation, label))


# def loadImages(folder, extension=None):
#     """
#     loads images into a list
#     - folder: all files in this folder will be added
#     - extension: specify a extension for the files.
#                  e.g. for '.jpg' use 'jpg'
#     """
#     imageList = []
#     filesInFolder = os.listdir(folder)
#     for f in filesInFolder:
#         if os.path.isfile(f):
#             if extension != None:
#                 if "." in f:
#                     if f.split(".")[-1].lower() == extension.lower():
#                         imageList.append(f)
#             else:
#                 imageList.append(f)
#     #fNameList = [f for f in os.listdir(capturePath) if os.path.isfile(os.path.join(capturePath, f)) and ".pcapng" in f]
#     return imageList
