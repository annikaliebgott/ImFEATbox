function Out = TopHatTrafoF(I,typeflag)
% Input:     - I: A 2D image
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.global: all features 
%                   + typeflag.form: all features
%                   + typeflag.transform: all features
%                   + typeflag.moments: only features based on moments
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x60) vector containing 60 metrics calculated from
%              the image transformed by top hat transformation
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

if ~exist('typeflag','var')
   typeflag.global = true; 
   typeflag.form = true;
   typeflag.transform = true;
   typeflag.moments = true;
end    


%% create morphological structruning element (STREL)
SE1 = strel('arbitrary', [1 0 0;1 0 0;1 0 1]);
SE2 = strel('ball', 15,5);
SE3 = strel('diamond', 13);
SE4 = strel('disk', 15);
SE5 = strel('line', 10, 45);
SE6 = strel('octagon', 3);
SE7 = strel('pair', [2 2]);
SE8 = strel('periodicline', 2, [1 -2]);
SE9 = strel('rectangle', [3 5]);
SE10 = strel('square', 11);

%% top-hat filtering
IM1 = double(imtophat(I,SE1));
IM2 = double(imtophat(I,SE2));
IM3 = double(imtophat(I,SE3));
IM4 = double(imtophat(I,SE4));
IM5 = double(imtophat(I,SE5));
IM6 = double(imtophat(I,SE6));
IM7 = double(imtophat(I,SE7));
IM8 = double(imtophat(I,SE8));
IM9 = double(imtophat(I,SE9));
IM10 = double(imtophat(I,SE10));

%% feature extraction

if (typeflag.global || typeflag.form || typeflag.transform)
% determine standard derivation
std_IM1 = std2(IM1);
std_IM2 = std2(IM2);
std_IM3 = std2(IM3);
std_IM4 = std2(IM4);
std_IM5 = std2(IM5);
std_IM6 = std2(IM6);
std_IM7 = std2(IM7);
std_IM8 = std2(IM8);
std_IM9 = std2(IM9);
std_IM10 = std2(IM10);

% calculate average value
mean_IM1 = mean2(IM1);
mean_IM2 = mean2(IM2);
mean_IM3 = mean2(IM3);
mean_IM4 = mean2(IM4);
mean_IM5 = mean2(IM5);
mean_IM6 = mean2(IM6);
mean_IM7 = mean2(IM7);
mean_IM8 = mean2(IM8);
mean_IM9 = mean2(IM9);
mean_IM10 = mean2(IM10);
end

% calculate 2nd and 4th moments
m_2_1 = moment(IM1,2);
m_2_2 = moment(IM2,2);
m_2_3 = moment(IM3,2);
m_2_4 = moment(IM4,2);
m_2_5 = moment(IM5,2);
m_2_6 = moment(IM6,2);
m_2_7 = moment(IM7,2);
m_2_8 = moment(IM8,2);
m_2_9 = moment(IM9,2);
m_2_10 = moment(IM10,2);

m_4_1 = moment(IM1,4);
m_4_2 = moment(IM2,4);
m_4_3 = moment(IM3,4);
m_4_4 = moment(IM4,4);
m_4_5 = moment(IM5,4);
m_4_6 = moment(IM6,4);
m_4_7 = moment(IM7,4);
m_4_8 = moment(IM8,4);
m_4_9 = moment(IM9,4);
m_4_10 = moment(IM10,4);

%average value of moments
mm_2_1 = mean2(m_2_1);
mm_2_2 = mean2(m_2_2);
mm_2_3 = mean2(m_2_3);
mm_2_4 = mean2(m_2_4);
mm_2_5 = mean2(m_2_5);
mm_2_6 = mean2(m_2_6);
mm_2_7 = mean2(m_2_7);
mm_2_8 = mean2(m_2_8);
mm_2_9 = mean2(m_2_9);
mm_2_10 = mean2(m_2_10);

mm_4_1 = mean2(m_4_1);
mm_4_2 = mean2(m_4_2);
mm_4_3 = mean2(m_4_3);
mm_4_4 = mean2(m_4_4);
mm_4_5 = mean2(m_4_5);
mm_4_6 = mean2(m_4_6);
mm_4_7 = mean2(m_4_7);
mm_4_8 = mean2(m_4_8);
mm_4_9 = mean2(m_4_9);
mm_4_10 = mean2(m_4_10);

% standard deviation of moments
ms_2_1 = std2(m_2_1);
ms_2_2 = std2(m_2_2);
ms_2_3 = std2(m_2_3);
ms_2_4 = std2(m_2_4);
ms_2_5 = std2(m_2_5);
ms_2_6 = std2(m_2_6);
ms_2_7 = std2(m_2_7);
ms_2_8 = std2(m_2_8);
ms_2_9 = std2(m_2_9);
ms_2_10 = std2(m_2_10);

ms_4_1 = std2(m_4_1);
ms_4_2 = std2(m_4_2);
ms_4_3 = std2(m_4_3);
ms_4_4 = std2(m_4_4);
ms_4_5 = std2(m_4_5);
ms_4_6 = std2(m_4_6);
ms_4_7 = std2(m_4_7);
ms_4_8 = std2(m_4_8);
ms_4_9 = std2(m_4_9);
ms_4_10 = std2(m_4_10);


%% return feature vector
if (typeflag.global || typeflag.form || typeflag.transform)
Out = [mean_IM1 mean_IM2 mean_IM3 mean_IM4 mean_IM5 mean_IM6 ...
    mean_IM7 mean_IM8 mean_IM9 mean_IM10 ...
    std_IM1 std_IM2 std_IM3 std_IM4 std_IM5 std_IM6 ...
    std_IM7 std_IM8 std_IM9 std_IM10 ...
    mm_2_1 mm_2_2 mm_2_3 mm_2_4 mm_2_5 mm_2_6 mm_2_7 mm_2_8 mm_2_9 mm_2_10 ...
    mm_4_1 mm_4_2 mm_4_3 mm_4_4 mm_4_5 mm_4_6 mm_4_7 mm_4_8 mm_4_9 mm_4_10 ...
    ms_4_1 ms_4_2 ms_4_3 ms_4_4 ms_4_5 ms_4_6 ms_4_7 ms_4_8 ms_4_9 ms_4_10 ...
    ms_2_1 ms_2_2 ms_2_3 ms_2_6 ms_2_7 ms_2_8 ms_2_9 ms_2_10 ...
    ms_2_5 ms_2_4
    ];
else
    Out = [mm_2_1 mm_2_2 mm_2_3 mm_2_4 mm_2_5 mm_2_6 mm_2_7 mm_2_8 mm_2_9 mm_2_10 ...
    mm_4_1 mm_4_2 mm_4_3 mm_4_4 mm_4_5 mm_4_6 mm_4_7 mm_4_8 mm_4_9 mm_4_10 ...
    ms_4_1 ms_4_2 ms_4_3 ms_4_4 ms_4_5 ms_4_6 ms_4_7 ms_4_8 ms_4_9 ms_4_10 ...
    ms_2_1 ms_2_2 ms_2_3 ms_2_6 ms_2_7 ms_2_8 ms_2_9 ms_2_10 ...
    ms_2_5 ms_2_4
    ];
end    
    

end