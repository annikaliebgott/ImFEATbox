function Out = HoughTrafoF(I,houghtype,arc_min,plotflag,typeflag)
% Input:     - I: A 2D image
%            - houghtype: a string defining which type of Hough transform
%              should be performed (circular, linear or both)
%              default: houghtype = 'both'
%            - arc_min: circles are treated relevant only if their arc is
%              at least arc_min. 
%              default: arc_min = pi/2
%            - plotflag: A locical flag to enable/disable visualization.
%              default: plotflag = false;
%            - typeflag: Struct of logicals to permit extracting features 
%              based on desired characteristics:
%                   + typeflag.global: all features 
%                   + typeflag.transform: all features 
%                   + typeflag.corr: only features based on correlation
%              default: all features are being extracted
%              For more information see README.txt
%
%
% Output:    - Out: A (1x393) vector containing 393 metrics calculated from
%              the Hough transform
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
   typeflag.transform = true;
   typeflag.form = true;
   typeflag.moments = true;
end    

if ~exist('houghtype','var')
   houghtype = 'both';
end

if ~exist('arc_min','var')
   arc_min = pi/2;
end    

if ~exist('plotflag','var')
   plotflag = false;
end    

if any(~isreal(I))
    I = real(I);
end

if strcmp(houghtype,'linear') || strcmp(houghtype,'both')
    
    %% Linear Hough transform
    % detect edges
    BW = edge(I, 'canny');     
    
    % linear Hough transform
    [H,theta,rho] = hough(BW);
    
    %% extract features
    length_theta = length(theta);
    theta_mean = mean(theta);
    theta_std = std(theta);
    length_rho = length(rho);
    rho_mean = mean(rho);
    rho_std = std(rho);
    H_std = std2(H);
    H_std2 = std(std(H));
    moment_2 = moment(H,2);
    moment_4 = moment(H,4);
    
    % find the peaks in the Hough transform matrix H
    P = houghpeaks(H, 5, 'threshold', ceil(0.3*max(H(:))));
    
    % find the lines in the image
    lines = houghlines(BW, theta, rho, P, 'FillGap', 5, 'MinLength', 7);
    N_lines = length(lines);
    
    % store starting and end points of all lines in array for further
    % calculations, calculate lengths of all lines
    line_start = zeros(N_lines,2);
    line_end = zeros(N_lines,2);
    length_lines = zeros(N_lines,1);
    for k = 1 : N_lines
        length_lines(k) = norm(lines(k).point1 - lines(k).point2);
        line_start(k,:) = lines(k).point1;
        line_end(k,:) = lines(k).point2;
    end
    
    % determine length of longest line
    if ~isempty(length_lines)
        max_length = max(length_lines);
    else
        max_length = 0;
    end    
    % calculate line with mean direction of all lines
    line_mean = [sum(line_start)./N_lines; sum(line_end)./N_lines];
    line_mean_slope = (line_mean(2,2)-line_mean(1,2))/(line_mean(2,1)-line_mean(1,1));
    
    % calculate length of line with mean direction of all lines
    line_mean_length = norm(line_mean(1,:) - line_mean(2,:));
    
    % calculat mean point of line with mean direction of all lines
    if N_lines > 1
        line_mean_point = mean(line_mean);
    else 
        line_mean_point = [0 0];
    end    
    
    % investigation of (normalized) parallelism of all detected lines by
    % counting the occurence of pairwise parallel lines (para_sum) and the
    % maximum amount of lines with the same gradient (para_max)
    para_sum = 0;
    grad = zeros(N_lines,1);
    for i = 1 : N_lines
        % calculate gradient of lines
        grad(i) = (line_end(i,2) - line_start(i,2))/(line_end(i,1) - line_start(i,1));
        for k = 1 : N_lines
            % count pairs of lines only once
            if ( i == k )
                break
            end
            % check for parallelism using the difference of the lines
            % gradients
            d_grad = abs(grad(k) - grad(i));
            if(d_grad < 0.01);
                para_sum = para_sum + 1;
            end
        end
    end
    
    % calculate the maximum amount of lines with the same gradient
    para_num = zeros(numel(unique(grad)),2);
    para_num(:,1) = unique(grad);
    for i = 1:length(grad)
        para_num(para_num(:,1) == grad(i),2) = para_num(para_num(:,1) == grad(i),2) + 1;
    end
    if ~isempty(para_num)
        para_max = max(para_num(:,2));
    else
        para_max = 0;
    end    
    
    % standard deviation of gradient
    line_grad_std = std(grad);
    
    % total length of all lines
    length_total = sum(length_lines);
    
    % mean length of all lines
    length_mean = mean(length_lines);
    
    % std of length of all lines
    length_std = std(length_lines);
    
    % total amount N_I of intersections of the lines in the image
    N_I = 0;
    for i = 1:N_lines
        for k = 1:N_lines
            % count pairs of lines only once
            if ( i == k )
                break
            end
            % check if lines are parallel
            if grad(i) ~= grad(k)
                Li = line_end(i,:) - line_start(i,:);
                Lk = line_end(k,:) - line_start(k,:);
                % calculate intersection of lines through the given points
                S = line_start(k,:)' + ((line_start(k,2)-line_start(i,2))*Li(1)/Li(2)...
                    - line_start(k,1)+line_start(i,1))/(Lk(1)-Lk(2)*Li(1)/Li(2))*Lk';
                
                % check whether intersection is within the lines start and
                % end points
                isIntersection = (((line_start(i,1) <= S(1) && S(1) <= line_end(i,1))...
                    || (line_start(i,1) >= S(1) && S(1) >= line_end(i,1)))...
                    && ((line_start(i,2) <= S(2) && S(2) <= line_end(i,2))...
                    ||(line_start(i,2) >= S(2) && S(2) >= line_end(i,2))))...
                    && (((line_start(k,1) <= S(1) && S(1) <= line_end(k,1))...
                    || (line_start(k,1) >= S(1) && S(1) >= line_end(k,1)))...
                    && ((line_start(k,2) <= S(2) && S(2) <= line_end(k,2))...
                    ||(line_start(k,2) >= S(2) && S(2) >= line_end(k,2)))) ;
                if isIntersection
                    N_I = N_I + 1;
                end
            end
            close all
        end
    end
    
    % average Number of intersections per line
    N_I_m = N_I/N_lines;
    
end

%% Circular Hough transform
if strcmp(houghtype,'circular') || strcmp(houghtype,'both')
    % Edge detection
    [~, threshold] = edge(I, 'sobel');
    fudgeFactor = 1;
    BWc = edge(I,'sobel', threshold * fudgeFactor);
    
    % Calculate circle centers
    centertable = zeros(0,5);
    radius = 10:100;
    circle = cell(length(radius),1);
    result = cell(length(radius),1);
    maxval_norm = zeros(length(radius),1);
    centers = cell(length(radius),1);
    result_norm = cell(length(radius),1);
    result_relevant = cell(length(radius),1);
    N_c = zeros(length(radius),1);
    N_c_relevant = zeros(length(radius),1);
    mean_arc = zeros(length(radius),1);
    mean_arc_relevant = zeros(length(radius),1);
    location = cell(length(radius),1);
    mean_loc = zeros(length(radius),2);
    std_loc = zeros(length(radius),2);
    
    
    for r = 1:length(radius)
        % define circle of radius radius(r) to filter the image with
        [x, y] = meshgrid(-radius(r):radius(r), -radius(r):radius(r));
        circle{r} = round(sqrt(x.^2 + y.^2)) == radius(r);
        
        % convolution of BWc with the circle defined in circle to get the
        % amount of non-zero pixels touched by circles of radius radius(r)
        % with mean (x,y) for each pixel (x,y) in BWc
        result{r} = conv2(double(BWc), double(circle{r}), 'same');
        % normalized results
        result_norm{r} = result{r}./radius(r);
        % save only relevant results representing a circle with an arc
        % greater than the minimal desired arc defined by arc_min
        result_relevant{r} = result_norm{r}(result_norm{r} >= arc_min);
        
        % mean location on x/y axis of the relevant circles with radius
        % radius(r)
        [x, y] = find(result_norm{r} >= arc_min);
        location{r} = [x y];
        mean_loc(r,:) = mean(location{r});
        std_loc(r,:) = std(location{r});
        
        % mean arc of detected circles with radius radius(r)
        mean_arc(r) = mean(result_norm{r}(:));
        mean_arc_relevant(r) = mean(result_relevant{r}(:));
        
        % Number of detected circles
        N_c(r) = nnz(result{r}(:));
        N_c_relevant(r) = nnz(result_relevant{r}(:));
        
        % calculate maximum number of non-zero pixels touched by circles of
        % radius radius(r)
        maxval_norm(r) = max(result_norm{r}(:));
        
    end
    
    %% extract features
    
    % total number of detected circles and number of relevant circles
    N_c_sum = sum(N_c);
    N_c_sum_relevant = sum(N_c_relevant);
    ratio_N_c = N_c_sum_relevant/N_c_sum;
    
    arc_mean = mean(mean_arc_relevant);
    
    % radius with the highest number of detected circles and relevant
    % circles. If more than one radius has the same highest number of
    % detected circles, the smallest radius is used as a feature.
    r_N_max = min(radius(N_c_relevant == max(N_c_relevant)));
    
    % mean radius for all detected circles and all relevant circles
    r_mean = radius*N_c_relevant/N_c_sum_relevant;
    
    % mean and std of maximum number of non-zero pixels per radius for all
    % detected circles and relevant circles
    maxval_mean = mean(maxval_norm);
    maxval_std = std(maxval_norm);
    
    % mean and std of location of relevant circles
    if N_c_sum_relevant > 1
        location_mean = mean(mean_loc);
        location_std = std(std_loc);
    else
        location_mean = [0 0];
        location_std = [0 0];
    end    
    
    %% visualization
    if plotflag
        for r = 1:length(radius)
            centers{r} = imregionalmax(result{r});
            [yc, xc] = find(centers{r} == 1);
            for k = 1:length(xc)
                centertable(end+1,:) = [result{r}(yc(k),xc(k)) xc(k) yc(k) r result{r}(yc(k),xc(k))/radius(r)];
            end
        end
        
        centertable = sortrows(centertable, -5);
        
        % chose only circles with circular arcs on detected edges that are
        % at least min_arc long
        centertable = centertable(centertable(:,5) >= arc_min,:);
        maxNrCircles = size(centertable,1);
        
        
        xc = centertable(:,2);
        yc = centertable(:,3);
        r = centertable(:,4);
        
%         figure; imshow(BWc);
        figure;
        I=mat2gray(I);
        imagesc(I)
        colormap(gray)
        hold on;
        plot(xc,yc,'x','LineWidth',2,'Color', 'yellow');
        for k = 1:maxNrCircles
            rectangle('Position',[xc(k)-radius(r(k)),yc(k)-radius(r(k)),2*radius(r(k)),2*radius(r(k))],'Curvature',[1,1],'EdgeColor', 'red', 'LineWidth',2)
        end
        hold off;
        
    end
end

%% return features
if strcmp(houghtype,'linear')
    if typeflag.global || typeflag.form || typeflag.transform      
        Out = [max_length line_mean_slope line_mean_length line_mean_point...
            para_sum para_max line_grad_std...
            length_total length_mean length_std...
            N_I N_I_m...
            length_theta theta_mean theta_std...
            length_rho rho_mean rho_std...
            H_std H_std2...
            moment_2 moment_4
            ];
    else
        Out = [moment_2 moment_4];
    end
end

if strcmp(houghtype, 'circular')
    Out = [N_c_sum N_c_sum_relevant ratio_N_c arc_mean r_N_max r_mean...
        maxval_mean maxval_std location_mean location_std];
end

if strcmp(houghtype, 'both')
    if typeflag.global || typeflag.form || typeflag.transform
        Out = [max_length line_mean_slope line_mean_length line_mean_point...
            para_sum para_max line_grad_std...
            length_total length_mean length_std...
            N_I N_I_m...
            length_theta theta_mean theta_std...
            length_rho rho_mean rho_std...
            H_std H_std2...
            moment_2 moment_4...
            N_c_sum N_c_sum_relevant ratio_N_c arc_mean r_N_max r_mean...
            maxval_mean maxval_std location_mean location_std];
    else
        Out = [moment_2 moment_4];
    end
end

