function Out = FormFactorF(I,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features
%              based on desired characteristics:
%                   + typeflag.global: all features
%                   + typeflag.form: all features
%                   + typeflag.corr: only features based on correlation
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x32) vector containing 32 metrics based on the form
%                   of detected objects in am image
%
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic
% and Interventional Radiology, University Hospital of Tuebingen, Germany
% and the Institute of Signal Processing and System Theory University of
% Stuttgart, Germany. Last modified: November 2016
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
% ************************************************************************


%% extract properties of regions in the image

if ~exist('typeflag','var')
    typeflag.global = true;
    typeflag.form = true;
    typeflag.corr = true;
end

% convert image
BW = im2bw(double(I));
[L,N] = bwlabel(BW);

% check if any objects have been detected
if N > 0
    % extract region properties
    s = regionprops(L, 'area', 'perimeter', ...
        'Orientation', 'Eccentricity', 'EquivDiameter', 'Solidity');
    
    % initialize variables
    roundness = zeros(1,N);
    form = zeros(1,N);
    area_obj = zeros(1,N);
    perimeter_obj = zeros(1,N);
    
    for i = 1 : N
        area_obj(i) = s(i).Area;
        perimeter_obj(i) = s(i).Perimeter;
        form(i) = 4*pi*area_obj(i)/(perimeter_obj(i)^2);
        roundness(i) = (4*s(i).Area )/(s(i).EquivDiameter*pi);
    end
    
    % avoid Inf values which cause problems in further calculations
    form(form == Inf) = max(form(form ~= Inf));
    
    eccentricity = [s.Eccentricity];
    orientation = [s.Orientation];
    solidity = [s.Solidity];
    
    
    
    %% feature extraction
    
    if (typeflag.global || typeflag.form)
        % mean
        mean_roundness = mean(roundness);
        mean_area = mean(area_obj);
        mean_perimeter = mean(perimeter_obj);
        mean_eccentricity = mean(eccentricity);
        mean_orientation = mean(orientation);
        mean_solidity = mean(solidity);
        mean_form = mean(form);
        
        % standard deviation
        std_roundness = std(roundness);
        std_area = std(area_obj);
        std_perimeter = std(perimeter_obj);
        std_eccentricity = std (eccentricity);
        std_orientation = std(orientation);
        std_solidity = std(solidity);
        std_form = std(form);
        
        % maximum/minimum values (orientation is neglected for this measure)
        max_roundness = max(roundness);
        max_area = max(area_obj);
        max_perimeter = max(perimeter_obj);
        max_eccentricity = max(eccentricity);
        max_solidity = max(solidity);
        max_form = max(form);
        
        min_roundness = min(roundness);
        min_area = min(area_obj);
        min_perimeter = min(perimeter_obj);
        min_eccentricity = min(eccentricity);
        min_solidity = min(solidity);
        min_form = min(form);
    end
    % calculate features from the correlation between form factor and roundness
    if N > 1
        corr_FR = xcorr(form, roundness, 'coeff');
        mean_corr = mean(corr_FR);
        std_corr = std(corr_FR);
        max_corr = max(corr_FR);
        min_corr = min(corr_FR);
        
        % determine cross correlation coefficient
        corrcoef_FR = corrcoef(form, roundness);
        corrcoef_FR = corrcoef_FR(1,2);
    else
        mean_corr = 0;
        std_corr = 0;
        max_corr = 0;
        min_corr = 0;
        corrcoef_FR = 0;
    end
    
    
    %% return feature vector
    if (typeflag.global || typeflag.form)
        Out = [N mean_eccentricity mean_solidity mean_orientation...
            mean_roundness mean_area mean_perimeter mean_form...
            std_eccentricity std_solidity std_orientation...
            std_roundness std_area std_perimeter std_form...
            max_eccentricity max_solidity max_roundness max_area max_perimeter max_form...
            min_eccentricity min_solidity min_roundness min_area min_perimeter min_form...
            corrcoef_FR mean_corr std_corr max_corr min_corr];
    else
        Out = [corrcoef_FR mean_corr std_corr max_corr min_corr];
    end
else
    Out = zeros(1,32);
end


end
