------------------------------------------------------------------------
1. How to use this toolbox
------------------------------------------------------------------------

There are two ways to use the functions provided by ImFEATbox:
1.) Directly use the feature extraction functions
    Note: in order to work properly, most of the feature extraction functions require preprocessing of the images,
    i.e. the images should be segmented and converted to a 255 gray scale. 
2.) Use the provided script 'ExtractFeatures.m' for feature extraction
    
------------------------------------------------------------------------
2. Usage of typeflag
------------------------------------------------------------------------

typeflag is a struct which serves to choose only certain categories of features for feature extraction. 
Note: the categories we assigned the features to are not necessarily the only ones the features might belong to. 
The following categories have been defined:

- typeflag.all: all features are being extracted
- typeflag.global: all global features are being extracted
- typeflag.local: all local features are being extracted
- typeflag.texture: all features depending on image texture are being extracted
- typeflag.transform: all features calculated from transformed images are being extracted
- typeflag.moments: all features calculated from moments are being extracted
- typeflag.corr: all features depending on correlation are being extracted
- typeflag.entropy: all features depending on entropy are being extracted

Most of the feature extraction algorithms provide measures fitting into more than one category. 
Also, there are algorithms where parts of the results belong to one category, other parts belong to others. 
As an example, the features calculated with the function IntensityF.m all belong to the categories 'global' 
and 'texture', but some measures additionally fit into the categories 'corr' and 'entropy'.

An up-to-date list of all features and the categories they have been assigned to can be found in features.pdf


************************************************************************
This file is part of ImFEATbox, a toolbox for image feature extraction and 
analysis. Implemented by the Department of Diagnostic and Interventional 
Radiology, University Hospital of Tuebingen, Germany and the Institute of 
Signal Processing and System Theory University of Stuttgart, Germany. 
Last modified: November 2016 

Available online at:
https://github.com/annikaliebgott/ImFEATbox

Contact: annika.liebgott@iss.uni-stuttgart.de
************************************************************************
