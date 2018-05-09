function Out = HarrisF(I,plotflag, N_s)
% Input:     - I: A 2D image
%            - plotflag: a flag to enable/disable visualization
%              Default: plotflag = false
%            - N_s: number of points considered to be the strongest points.
%              Default: N_s = ceil(0.1*N) (N: number of detected corners)
%
% Output:    - Out: A (1x10) vector containing 10 metrics calculated based
%                   on corner points detected with Harris-Stephens algorithm
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
if~exist('plotflag','var')
    plotflag = false;
end

% convert image
I = double(I);

if any(~real(I))
    I = real(I);
end

% calculation of corner points using Harris-Stephens algorithm
corners = detectHarrisFeatures(I);


%% extract features
c = double(corners.Location);
N = double(corners.Count);

if N > 0
    
    if ~exist('N_s','var') || N_s > 0.5*N
        N_s = ceil(0.2*N);
    end
    
    c_s = double(corners.selectStrongest(N_s).Location);
    
    % determine center  of gravity for all points
    x_gravity = sum(c(:,1)./N);
    y_gravity = sum(c(:,2)./N);
    
    % determine center of gravity for the strongest points
    x_gravity_s = sum(c_s(:,1)./N_s);
    y_gravity_s = sum(c_s(:,2)./N_s);
    
    % display the results
    if plotflag
        %     imshow(I);
        imagesc(I);
        colormap(gray);
        hold on;
        scatter(c(:,1),c(:,2),'+');
        scatter(c_s(:,1),c_s(:,2),'g+');
        plot(x_gravity_s, y_gravity_s,'r*','MarkerSize',10);
    end
    
    % density of corner points
    dens = N / (numel(I));
    
    %standard derivation
    std_x = std(c(:,1));
    std_y = std(c(:,2));
    std_x_s = std(c_s(:,1));
    std_y_s = std(c_s(:,2));
    
    
    %% return feature vector
    Out = [N x_gravity y_gravity x_gravity_s y_gravity_s dens...
        std_x std_y std_x_s std_y_s];
else
    Out =  zeros(1,10);
end
end
