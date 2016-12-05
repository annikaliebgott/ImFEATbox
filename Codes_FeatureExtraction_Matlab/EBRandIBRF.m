    function Out = EBRandIBRF(I,typeflag)
    % Input:     - I: A 2D image
    %            - typeflag: Struct of logicals to permit extracting features 
    %              based on desired characteristics:
    %                   + typeflag.local: all features 
    %                   + typeflag.form: all features
    %                   + typeflag.texture: all features
    %                   + typeflag.corr: only features based on correlation
    %              default: all features are being extracted
    %              For more information see README.txt
    %
    %
    % Output:    - Out: A (1x1211) vector containing 1211 metrics 
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
    % Implementation based on:  Tinne Tuytelaars and Luc Van Gool. 2004. 
    %                           Matching Widely Separated Views Based on Affine 
    %                           Invariant Regions. Int. J. Comput. Vision 59, 1 
    %                           (August 2004), 61-85.


    if ~exist('typeflag','var')
       typeflag.local = true; 
       typeflag.form = true;
       typeflag.texture = true;
       typeflag.corr = true;
    end    

    N2 = size(I,2);
    N1 = size(I,1);

    % convert image to a even sized image
    if (mod(N1,2) == 0 && mod(N2,2)==0 )

    elseif(mod(N1,2)>0)
        I = I(1:N1-1,:);
    elseif(mod(N2,2)>0)
        I = I(:,1:N2-1);
    end


    %% Detect Corners and Edges
    % select Harris corners
    C1 = corner(I);

    % find Canny edges
    BW = edge(double(I),'Canny');

    % check if corners lay on edges
    temp = [];
    for i = 1 : size(C1,1)
        if C1(i,1) < size(BW,1) && C1(i,2) < size(BW,2)
            if (BW(C1(i,1), C1(i,2)) == 1)
                temp = [temp; [C1(i,1) C1(i,2)]];
            end
        end
    end


    % Number of corners
    N_corners = size(temp,1);


    %% construct the smallest parallelogram with the points

    % matrix with all paralleograms
    para = [];
    para_total=[];
    orientation_para = [];
    area_total = 0;
    perim_total = 0;

    % loop over all coners in temp
    for i = 1 : N_corners
        % extract region of interesst through the corner point
        % Specify the locations of objects in the image using row and column indices.
        c = temp(i,1);
        r = temp(i,2);
        BW2 = bwselect(BW,c,r,8);

        % determine indices of non-zero elements
        [y, x] = find(BW2);

        if (length(x) > 2)

            % determine the minimum bounding parallogram
            [pgx,pgy,area,perim] = minboundparallelogram(x, y);

            % for testing
            PlotFlag = 0;
            if PlotFlag
                imshow(I), hold on;  plot(pgx(:),pgy(:),'LineWidth',2,'Color','green');
            end

            para = [para; [pgx,pgy]];
            para_total = [para_total; para]; %new: ;

            area_total = area_total + (area/N_corners);
            perim_total = perim_total + (perim/N_corners);

            % determine the slope and degree ot the slope
            m_slope = abs((pgy(1,1) - pgy(2,1)) / (pgx(1,1) - pgx(2,1)));
            orientation_para = [orientation_para tan(m_slope)];

        end
    end


    %% extract features

    % number of paralleograms
    num_para = size(para,1);

    if size(para_total,1) < 1 && size(para_total,2)
        % to maintain to translation invariant feature
        para_total = para_total/(size(I,1)+size(I,2));
        para_x = shiftdim(para_total(:,1),1);
        para_y = shiftdim(para_total(:,2),1);

        % determine the orientation of the parallelogramms
        orientation_para = anglemean(orientation_para);
    else
        para_x = 0;
        para_y = 0;
        orientation_para = 0;
    end

    C3 = cornermetric(I, 'MinimumEigenvalue');

    % use different sensitivity factors
    z = (0.01:0.01:0.2);
    % initiate variables
    u_tr = zeros(size(z));
    u_rank = zeros(size(z));
    u_det = zeros(size(z));
    coefs = zeros(1,4*size(z,2));
    u_mean = zeros(size(z));
    u_std = zeros(1,2*size(z,2));
    corr_C2C3 =  zeros(1,size(z,2));
    eU = zeros(1,50*size(z,2));

    for i = 1:size(z,2)
        C2 = cornermetric(I, 'Harris', 'SensitivityFactor', z(i));

        if (typeflag.local || typeflag.form || typeflag.texture)
            % use the cornermetic matrices C2 and C3
            U = C3*transpose(C2);

            u_tr(i) = trace(U);
            u_rank(i) = rank(U);
            u_det(i) = det(U);
            u_mean(i) = mean2(U);
            u_std(2*(i-1)+1:2*i) = [std2(U) std(std(U))];
            help = shiftdim(eig(U),1);
            eU(50*(i-1)+1:50*i) = help(1:50);


            if (size(U,1) > size(U,2))
                maxSize = size(U,2);
            else
                maxSize = size(U,1);
            end


        % use 4 fixed coefficients as features
        coefs(4*(i-1)+1:4*i) = [U(round(maxSize/2),round(maxSize/2))...
            U(round(maxSize/4),round(maxSize/4)) U(round(maxSize/2),round(maxSize/4))...
            U(round(maxSize/4),round(maxSize/2))];
        end

        % correlation
        Corr = corrcoef(C2,C3);
        corr_C2C3(i) = Corr(1,2);
    end



    %% return feature vector

    if ~(typeflag.local || typeflag.form || typeflag.texture)
        Out = [mean(corr_C2C3) max(corr_C2C3) min(corr_C2C3) std(corr_C2C3)];
    else
        Out = [N_corners u_tr u_rank u_det u_mean u_std...
            mean(corr_C2C3) max(corr_C2C3) min(corr_C2C3) std(corr_C2C3)...
            eU coefs area_total perim_total para_x para_y num_para orientation_para];
    end
    end

