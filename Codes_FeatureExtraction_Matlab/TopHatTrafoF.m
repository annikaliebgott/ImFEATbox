function Out = TopHatTrafoF(I,SE,typeflag)
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

if ~exist('SE','var')
    % create morphological structruning element (STREL)
    SE = cell(zeros);
    SE{1} = strel('arbitrary', [1 0 0;1 0 0;1 0 1]);
    SE{2} = strel('ball', 15,5);
    SE{3} = strel('diamond', 13);
    SE{4} = strel('disk', 15);
    SE{5} = strel('line', 10, 45);
    SE{6} = strel('octagon', 3);
    SE{7} = strel('pair', [2 2]);
    SE{8} = strel('periodicline', 2, [1 -2]);
    SE{9} = strel('rectangle', [3 5]);
    SE{10} = strel('square', 11);
end


%% top-hat filtering
N_SE = size(SE,2);
IM = cell(zeros);
for i = 1:N_SE
    IM{i} = double(imtophat(I,SE{i}));
    % remove possible Inf values
    IM{i}(IM{i} == Inf) = max(IM{i}(IM{i}(:)<Inf));
end

%% feature extraction

if (typeflag.global || typeflag.form || typeflag.transform)
    % determine standard derivation and average value
    std_IM = zeros(1,N_SE);
    mean_IM = zeros(1,N_SE);
    for i = 1:N_SE
        std_IM(i) = std2(IM{i});
        mean_IM(i) = mean2(IM{i});
    end
end

% calculate 2nd and 4th moments
m_2 = zeros(N_SE,size(I,2));
    m_4 = zeros(N_SE,size(I,2));
    mm_2 = zeros(1,N_SE);
    mm_4 = zeros(1,N_SE);
    ms_2 = zeros(1,N_SE);
    ms_4 = zeros(1,N_SE);

for i = 1:N_SE
    m_2(i,:) = moment(IM{i},2);
    m_4(i,:) = moment(IM{i},4);

    %average value of moments
    mm_2(i) = mean(m_2(i,:));   
    mm_4(i) = mean(m_4(i,:));
  
    % standard deviation of moments
    ms_2(i) = std(m_2(i,:));
    ms_4(i) = std(m_4(i,:));

end

%% return feature vector
if (typeflag.global || typeflag.form || typeflag.transform)
    Out = [mean_IM std_IM mm_2 mm_4 ms_2 ms_4];
else
    Out = [mm_2 mm_4 ms_2 ms_4];
end


end