To start feature reduction, adapt and run the script startFeatureReduction.m 

In case of unsupervised reduction algorithms, a feature matrix of size [N_samples x N_features] must be provided. This is the same format as feat_vector, the output of the feature extraction framework of this toolbox. To use the supervised reduction algorithms, a label vector of size [N_samples x 1] corresponding to the feature matrix must be provided additionally.

Part of the implemented algorithms are based on or extracted from the Matlab Toolbox for Dimensionality Reduction by Laurens van der Maaten (http://homepage.tudelft.nl/19j49)
(c) Laurens van der Maaten, Delft University of Technology

Part of the implemented algorithms are based on or extracted from the MatlabFunc toolbox by Deng Cai 
(c) Deng Cai, Zhejiang University
