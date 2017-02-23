------------------------------------------------------------------------
1. How to use this toolbox
------------------------------------------------------------------------

There are two ways to use the functions provided by ImFEATbox:
1.) Directly use the feature extraction functions
    Note: in order to work properly, most of the feature extraction functions require preprocessing of the images,
    i.e. the images need to be segmented and/or converted to N gray scales (done inside the script). 
2.) Use the provided function 'fExtractFeatures.m' for feature extraction (see example in example_extractFeatures.m)
    
------------------------------------------------------------------------
2. Feature categories
------------------------------------------------------------------------

Note: 
The categories we assigned the features to are not necessarily the only ones the features might belong to. 
Most of the feature extraction algorithms provide measures fitting into more than one category. 
Also, there are algorithms where parts of the results belong to one category, other parts belong to others. 
As an example, the features calculated with the function IntensityF.m all belong to the categories 'global' 
and 'texture', but some measures additionally fit into the categories 'corr' and 'entropy'.

- 'all': all features are being extracted
- 'global': all global features are being extracted
- 'local': all local features are being extracted
- 'texture': all features depending on image texture are being extracted
- 'transform': all features calculated from transformed images are being extracted
- 'moments': all features calculated from moments are being extracted
- 'corr': all features depending on correlation are being extracted
- 'entropy': all features depending on entropy are being extracted


An up-to-date list of all features and the categories they have been assigned to can be found in features.pdf
OR from the function itself:
cFeatures = fExtractFeatures([],[], '-getFeatures');


************************************************************************
This file is part of ImFEATbox, a toolbox for image feature extraction and 
analysis. Implemented by the Department of Diagnostic and Interventional 
Radiology, University Hospital of Tuebingen, Germany and the Institute of 
Signal Processing and System Theory University of Stuttgart, Germany. 
Last modified: February 2017

Available online at:
https://github.com/annikaliebgott/ImFEATbox

Contact: annika.liebgott@iss.uni-stuttgart.de
		 thomas.kuestner@iss.uni-stuttgart.de
************************************************************************