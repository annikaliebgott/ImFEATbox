function [hFeatures, cUsedFeatures, typeflag] = fParseFeatureAlgos( cFeatureAlgosIn, sFeaturePara )
% prepare feature algorithms which should be run
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic
% and Interventional Radiology, University Hospital of Tuebingen, Germany
% and the Institute of Signal Processing and System Theory University of
% Stuttgart, Germany. Last modified: February 2017
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de, thomas.kuestner@iss.uni-stuttgart.de
% ************************************************************************

%% get feature parameters
% -> to retrieve correct size of some features
currpath = fileparts(mfilename('fullpath'));
[pathPara, filenamePara, ~] = fileparts(sFeaturePara);
cd(pathPara);
eval([filenamePara,';']);
cd(currpath);

%% implemented features
% ATTENTION:
% nFeatures is determined by set parametrization (for some) -> future: invoke call to function handle with parametrization to return size of feature vector for current parametrization
% grouptypes: binary encoded with (undefined,global,local,corr,gradientGroup,moments,texture,form,entropy,transform)
%                    algorithm name      function handle     input parameter     nFeatures  group   subgroup   grouptypes
cMapFeatureGroup = { % GLOBAL FEATURES                                                      0                  (undefined,global, local,corr, gradientGroup,moments, texture,form, entropy,transform)
                     % #################################################################
                     % -----------------------------------------------------------------
                     % Intensity features                                                           0
                     % -----------------------------------------------------------------
                     'intensity',        'IntensityF',       {'typeflag'},          7,      0,      0,         bin2dec('01 01 00 10 10');
                     'histogram',        'HistogramF',       {'typeflag'},          6,      0,      0,         bin2dec('01 00 00 10 10');
                     'svd',              'SVDF',             {},                    780,    0,      0,         bin2dec('01 00 00 10 00');
                     'gradient',         'GradientF',        {'typeflag','gradtype'},81,    0,      0,         bin2dec('01 00 10 10 10');
                     % -----------------------------------------------------------------
                     % Geometrical features                                                         1
                     % -----------------------------------------------------------------
                     'glcm',             'GLCMF',            {'GLCMParameters', 'typeflag'}, 672, 0, 1,        bin2dec('01 01 00 11 10');
                     'runlength',        'RunLengthF',       {},                     44,     0,      1,        bin2dec('01 00 00 10 00');
                     'fractaldimension', 'FractalDimensionF',{'plotflag','width'},   (log(width)/log(2))*3+3,     0,      1,        bin2dec('01 00 00 01 00');
                     'formfactor',       'FormFactorF',      {'typeflag'},           32,     0,      1,        bin2dec('01 01 00 01 00');
                     % -----------------------------------------------------------------
                     % Transformation features                                              2
                     % -----------------------------------------------------------------
                     'fourier',          'FourierTrafoF',   {'typeflag'},            300,     0,       2,      bin2dec('01 01 01 00 01');
                     'dct',              'DCTF',            {'typeflag'},            2901,    0,       2,      bin2dec('01 01 00 00 01');
                     'hankel',           'DHankelF',        {'typeflag'},            75,      0,       2,      bin2dec('01 01 00 00 01');
                     'distancetrafo',    'DistanceTrafoF',  {'typeflag'},            56,      0,       2,      bin2dec('01 01 00 00 01');
                     'tophat',           'TopHatTrafoF',    {'SE','typeflag'},       6*size(SE,2), 0,       2, bin2dec('01 00 01 01 01');
                     'skeletonization',  'SkeletonizationF',{'typeflag'},            17,      0,       2,      bin2dec('01 01 00 00 01');
                     'unitarytrafo',     'UnitaryTrafoF',   {'transformation'},      73,      0,       2,      bin2dec('01 00 00 00 01');
                     'hough',            'HoughTrafoF',     {'houghtype','arc_min','plotflag','typeflag'}, 393, 0, 2, bin2dec('01 00 01 01 01');
                     'gaborfilter',      'GaborFilterF',    {'typeflag','gradtype','scale','orientation','plotflag'},3600, 0,2, bin2dec('01 00 10 10 11');
                     'wavelet',          'WaveletTrafoF',   {'typeflag','wavelettype'},4759,    0,     2,        bin2dec('01 01 01 00 11');
                    % -----------------------------------------------------------------
                    % Moment features                                                       3
                    % -----------------------------------------------------------------
                     'zernike',          'ZernikeF',        {},                       92,     0,      3,       bin2dec('01 00 01 00 00');
                     'hu',               'HuF',             {},                       8,      0,      3,       bin2dec('01 00 01 00 00');
                     'affinemoments',    'AffineMomentsF',  {},                       6,      0,      3,       bin2dec('01 00 01 00 00');

                    % LOCAL FEATURES                                                  1
                    % #################################################################
                    % -----------------------------------------------------------------
                    % Region features                                                         0
                    % -----------------------------------------------------------------
                     'lbp',              'LocalBinaryPatternF', {'RadialLBP'},        1024,   1,      0,       bin2dec('00 10 00 10 00');
                     'mser',             'MSERF',           {'plotflag'},             15,      1,      0,      bin2dec('00 10 00 10 00');
                     'saliency',         'SalientRegionF',  {'border_Salient','R1','typeflag'},8,1, 0,         bin2dec('00 10 01 10 00');
                     'quadtree',         'QuadtreeDecompositionF', {'threshold_quadtree'},5, 1, 0,             bin2dec('00 10 00 10 00');
                     'ebribr',           'EBRandIBRF',      {'typeflag'},             1211,   1,      0,       bin2dec('00 11 00 11 00');
                     'connectivity',     'ConnectivityF',   {},                       47,     1,      0,       bin2dec('00 10 00 10 00');
                     'sector',           'SectorF',         {},                       5,     1,      0,        bin2dec('00 10 00 00 00');
                     'lacunarity',       'LacunarityF',     {},                       6,    1,      0,         bin2dec('00 11 00 10 00');
                    % -----------------------------------------------------------------
                    % Line features                                                           1
                    % -----------------------------------------------------------------
                     'harris',           'HarrisF',         {'plotflag','N_s_Harris'},10,     1,      1,       bin2dec('00 10 00 10 00');
                     'line',             'LineProfileF',    {'plotflag','typeflag'},  122,    1,      1,       bin2dec('00 11 01 10 00');
                    % -----------------------------------------------------------------
                    % Point features                                                          2
                    % -----------------------------------------------------------------
                     'law',              'LawF',            {},                       58,     1,      2,       bin2dec('00 10 01 10 00');
                     'log',              'LoGF',            {'sigma','N_blobs','typeflag'},261,   1, 2,        bin2dec('00 10 01 10 00');
                     'gilles',           'GillesF',         {'radius_Gilles','threshold_Gilles'},6,1,2,        bin2dec('00 10 00 00 10');

                    % FEATURE DESCRIPTORS                                             2
                    % #################################################################
                     'surf',             'SURF',            {'N_s_SURF'},             11,    2,     0,         bin2dec('00 00 00 10 00');
                     'losib',            'LOSIBF',          {'radius','neighbors'},   34,    2,     0,         bin2dec('00 00 00 10 00');
                     'rcovd',            'RCovDsF',         {'typeflag'},             15,    2,     0,         bin2dec('00 01 01 00 00')};

% adapt size of some parameters
if transformation.WH && transformation.H && transformation.R && transformation.C
    cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 73;
else
    if transformation.WH
        if transformation.H
            if transformation.R
                % WH + H + R
                cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 15;
            elseif transformation.C
                % WH + H + C
                cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 68;
            else
                % WH + H
                cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 10;
            end
        elseif transformation.R
            if transformation.C
                % WH + R + C
                cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 68;
            else
                % WH + R
                cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 10;
            end
        elseif transformation.C
            % WH + C
            cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 63;
        else
            % WH
            cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 5;
        end
    elseif transformation.H
        if transformation.R
            if transformation.C
                % H + R + C
                cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 68;
            else
                % H + R
                cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 10;
            end
        elseif transformation.C
            % H + C
            cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 63;
        else
            % H
            cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 5;
        end
    elseif transformation.R
        if transformation.C
            % R + C
            cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 63;
        else
            % R
            cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 5;
        end
    else
        % C
        cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'unitarytrafo'),4} = 58;
    end
end
if strcmp(houghtype,'both')
    cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'hough'),4} = 393;
end
if strcmp(houghtype,'linear')
    cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'hough'),4} = 381;
end
if strcmp(houghtype,'circular')
    cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),'hough'),4} = 12;
end


%% parse inputs
cGroups = cMapFeatureGroup(:,end);
cFeatureAlgoGroups = [];
cFeatureAlgos = cell(size(cFeatureAlgosIn));
for iI = 1:length(cFeatureAlgosIn)
    switch lower(cFeatureAlgosIn{iI})
        case 'all'
            cFeatureAlgos = cMapFeatureGroup(:,1);
            break;
        case 'global'
            lMask = cellfun(@(x) bitand(x,2^8) > 0, cGroups);
            cFeatureAlgoGroups = cat(1, cFeatureAlgoGroups, cMapFeatureGroup(lMask,1));
        case 'local'
            lMask = cellfun(@(x) bitand(x,2^7) > 0, cGroups);
            cFeatureAlgoGroups = cat(1, cFeatureAlgoGroups, cMapFeatureGroup(lMask,1));
        case 'corr'
            lMask = cellfun(@(x) bitand(x,2^6) > 0, cGroups);
            cFeatureAlgoGroups = cat(1, cFeatureAlgoGroups, cMapFeatureGroup(lMask,1));
        case 'gradientGroup'
            lMask = cellfun(@(x) bitand(x,2^5) > 0, cGroups);
            cFeatureAlgoGroups = cat(1, cFeatureAlgoGroups, cMapFeatureGroup(lMask,1));
        case 'moments'
            lMask = cellfun(@(x) bitand(x,2^4) > 0, cGroups);
            cFeatureAlgoGroups = cat(1, cFeatureAlgoGroups, cMapFeatureGroup(lMask,1));
        case 'texture'
            lMask = cellfun(@(x) bitand(x,2^3) > 0, cGroups);
            cFeatureAlgoGroups = cat(1, cFeatureAlgoGroups, cMapFeatureGroup(lMask,1));
        case 'form'
            lMask = cellfun(@(x) bitand(x,2^2) > 0, cGroups);
            cFeatureAlgoGroups = cat(1, cFeatureAlgoGroups, cMapFeatureGroup(lMask,1));
        case 'entropy'
            lMask = cellfun(@(x) bitand(x,2^1) > 0, cGroups);
            cFeatureAlgoGroups = cat(1, cFeatureAlgoGroups, cMapFeatureGroup(lMask,1));
        case 'transform'
            lMask = cellfun(@(x) bitand(x,2^0) > 0, cGroups);
            cFeatureAlgoGroups = cat(1, cFeatureAlgoGroups, cMapFeatureGroup(lMask,1));
        case '-getfeatures'
            hFeatures = []; typeflag = [];
            cUsedFeatures = cMapFeatureGroup(:,[1 4]);
            return;
        otherwise
            if(~any(strcmp(cMapFeatureGroup(:,1),cFeatureAlgosIn{iI})))
                error('Selected feature is not available');
            end
            cFeatureAlgos{iI} = cFeatureAlgosIn{iI};
    end
end
if(~isempty(cFeatureAlgos{1}))
	cFeatureAlgos = cat(1,cFeatureAlgos,cFeatureAlgoGroups);
else
	cFeatureAlgos = cFeatureAlgoGroups;
end

%% set typeflags
bGroups = 0;
for iI = 1:length(cFeatureAlgos)
    bGroups = bitor(bGroups,cMapFeatureGroup{strcmp(cMapFeatureGroup(:,1),cFeatureAlgos{iI}),7});
end
typeflag = struct;
typeflag.global = false;
typeflag.local = false;
typeflag.corr = false;
typeflag.gradient = false;
typeflag.moments = false;
typeflag.texture = false;
typeflag.form = false;
typeflag.entropy = false;
typeflag.transform = false;
for iI=0:8
    if(bitand(bGroups,2^iI) > 0)
        switch iI
            case 8
                typeflag.global = true;
            case 7
                typeflag.local = true;
            case 6
                typeflag.corr = true;
            case 5
                typeflag.gradient = true;
            case 4
                typeflag.moments = true;
            case 3
                typeflag.texture = true;
            case 2
                typeflag.form = true;
            case 1
                typeflag.entropy = true;
            case 0
                typeflag.transform = true;
        end
    end
end


%% get handles
hFeatures = cell(length(cFeatureAlgos),1);
cUsedFeatures = cell(length(cFeatureAlgos),size(cMapFeatureGroup,2));
for iI = 1:length(cFeatureAlgos)
    iInd = find(strcmp(cMapFeatureGroup(:,1),cFeatureAlgos{iI}));
    cUsedFeatures(iI,:) = cMapFeatureGroup(iInd,:);
    if(isempty(cMapFeatureGroup{iInd,3}))
        eval(sprintf('hFeatures{iI} = @(x) %s(x);', cMapFeatureGroup{iInd,2}));
    else
        sPara = '';
        for iPara=1:length(cMapFeatureGroup{iInd,3})
            sPara = [sPara, cMapFeatureGroup{iInd,3}{iPara}, ','];
        end
        sPara = sPara(1:end-1);
        eval(sprintf('hFeatures{iI} = @(x) %s(x,%s);', cMapFeatureGroup{iInd,2},sPara));
    end
end



end
