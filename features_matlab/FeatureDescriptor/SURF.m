function Out = SURF(I,N_s,plotflag)
% Input:     - I: A 2D image
%            - N_s: desired number of strongest points to be selected.
%              Default: N_s = 25
%
%
% Output:    - Out: A (1x11) vector containing 11 metrics
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
%
% Implementation based on:  Herbert Bay, Andreas Ess, Tinne Tuytelaars,
%                           and Luc Van Gool. 2008. Speeded-Up Robust
%                           Features (SURF). Comput. Vis. Image Underst.
%                           110, 3 (June 2008), 346-359


%% calculation of SURF key points

% convert image
I = double(I);

if any(~real(I))
    I = real(I);
end

% calculation of the key points
points = detectSURFFeatures(I);

%Number of detected points
N_points = double(points.length);

if N_points > 0
    % Location of all key points
    L = double(points.Location);

    if ~exist('N_s','var') || N_s > 0.5*N_points
        N_s = ceil(0.2*N_points);
    end
    
    
    % Location of the N_s strongest key points
    L_s = double(points.selectStrongest(N_s).Location);
    
    % metric describing the strength of the detected points
    M = double(points.Metric);
    
    
    
    %% feature extraction
    
    % mean of M
    mean_M = mean(M);
    
    % determine center  of gravity for all points
    x_gravity = sum(L(:,1))./N_points;
    y_gravity = sum(L(:,2))./N_points;
    
    % determine center of gravity for the strongest points
    x_gravity_s = sum(L_s(:,1))./N_s;
    y_gravity_s = sum(L_s(:,2))./N_s;
    
    % calculate standard deviation of M, L and L_s
    std_M = std(M);
    
    if N_points > 1
        std_L = std(L);
    else
        std_L = [0 0];
    end
    
    if N_s > 1
        std_L_s = std(L_s);
    else
        std_L_s = [0 0];
    end
    
    %% plot selected points if plotflag is set to 1
    if plotflag
        imshow(I); hold on;
        plot(points.selectStrongest(N_s)); hold all;
        plot(x_gravity_s,y_gravity_s,'*r'); hold all;
        plot(x_gravity, y_gravity,'*r');
    end
    
else
    std_M = 0;
    std_L = [0 0];
    std_L_s = [0 0];
    mean_M = 0;
    x_gravity = 0;
    y_gravity = 0;
    x_gravity_s = 0;
    y_gravity_s = 0;
end

%% return feature vector
Out = [N_points std_M std_L std_L_s...
    mean_M x_gravity y_gravity x_gravity_s y_gravity_s];

