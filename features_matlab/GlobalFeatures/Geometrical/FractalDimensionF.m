function Out = FractalDimensionF(I,plotflag,width) 
% Input:     - I: A 2D image
%            - plotflag: A locical flag to enable/disable visualization.
%              Default: false
%            - width: largest size of the box. Default: width = 256
%
% Output:    - Out: A (1 x (log(width)/log(2))*3 + 3) vector containing
%                   metrics based on self-similarity of image structures
%
%
% Implemented algorithms:   1. Box-Counting(Haus) (BC)
%                           2. Differential Box-Counting (MBC)
%                           3. Triangular Prism Surface Area (TPSA)
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

if ~exist('plotflag','var')
    plotflag = false; 
end    
if ~exist('width','var')
    width = 256;
end

dim = ndims(I); 
if dim > 2
    error('Maximum dimension should be 2.');
end


% Spliting the Image into the small grids and copute the fractal dimension

% The largest size of the box
p = log(width)/log(2);   
RescaledI = imresize(I,[width width]);

% Allocation of the number of box size
counter=0;
counter_dbc = 0;
counter_tpsa = 0;
step = width./2.^(1:p);
N_mbc = zeros(1,length(step));
N_tpsa = zeros(1,length(step));
N_b = zeros(1,length(step));

%% 2D boxcount 
for n = 1:length(step)
    stepnum = step(n);
    for i = 1:stepnum:width 
        for j = 1:stepnum:width
            
            % Get the Grid in each level
            testim = RescaledI(i:i+stepnum-1,j:j+stepnum-1);
            
            % Differential(Modified) Box Counting
            MaxGrayLevel = max(testim(:));
            MinGrayLevel = min(testim(:));
            GridCont = MaxGrayLevel-MinGrayLevel+1;
            counter_dbc = counter_dbc + GridCont;
            % Differential(Modified) Box Counting (MBC)
            
            %Triangular Prism Surface Area (TPSA)
            a = testim(1,1);
            b = testim(1,end);
            c = testim(end,1);
            d = testim(end,end);
            e = (a+b+c+d)/4;
            
            w = sqrt(((b-a)^2) + (stepnum^2));
            x = sqrt(((c-b)^2) + (stepnum^2));
            y = sqrt(((d-c)^2) + (stepnum^2));
            z = sqrt(((a-d)^2) + (stepnum^2));
            
            o = sqrt(((a-e)^2) + (0.5*stepnum^2));
            p2 = sqrt(((b-e)^2) + (0.5*stepnum^2));
            q = sqrt(((c-e)^2) + (0.5*stepnum^2));
            t = sqrt(((d-e)^2) + (0.5*stepnum^2));
            
            % Using Herons Formula
            sa = (w+p2+o)/2;
            sb = (x+p2+q)/2;
            sc = (y+q+t)/2;
            sd = (z+o+t)/2;
            
            % Areas of Traiangle
            S_ABE = sqrt(sa*(sa-w)*(sa-p2)*(sa-o));
            S_BCE = sqrt(sb*(sb-x)*(sb-p2)*(sb-q));
            S_CDE = sqrt(sc*(sc-q)*(sc-t)*(sc-y));
            S_DAE = sqrt(sd*(sd-z)*(sd-o)*(sd-t));
            SurfaceArea = S_ABE + S_BCE + S_CDE + S_DAE;
            counter_tpsa = counter_tpsa + SurfaceArea;
            %Triangular Prism Surface Area

            
            % Basic Box Counting (BC)
            if (size(find(testim~=0),1)~=0)
                counter = counter+1;
            end	
            
        end
    end
    N_mbc (1,n) = counter_dbc;
    N_tpsa (1,n) = counter_tpsa;
    N_b(1,n) = counter;
    counter = 0;
    counter_dbc = 0;
    counter_tpsa = 0;
end

% Resolution
r0 = 2.^(p:-1:1);

% Dimension of BC
x0 = log(r0);
y0 = log(N_b);
FDMat_BC = y0./x0;
D0 = polyfit(x0, y0, 1);
FD_BC = D0(1);

% Dimension of MBC
x1 = log(r0);
y1 = log(N_mbc);
FDMat_MBC = y1./x1;
D1 = polyfit(x1, y1, 1);
FD_MBC = D1(1);

% Dimension of TPSA
x2 = log(r0);
y2 = log(N_tpsa);
FDMat_TPSA = y2./x2;
D2 = polyfit(x2, y2, 1);
FD_TPSA = 2 - D2(1);

% Plotting
if plotflag
    
    % Figure 1
    f0 = polyval(D0,x0);
    figure, plot(x0,y0,'-*','color','b','LineWidth',1.5)
    grid
    hold on, plot(x0,f0,'-*','color','k','LineWidth',1.5)
    legend('The FD Line','The Best Fitted Line','Location','NorthEast')
    xlabel('log(r)')
    ylabel('log(N)')
    
    % Figure 2
    f1 = polyval(D1,x1);
    figure, plot(x1,y1,'-*','color','b','LineWidth',1.5)
    grid
    hold on, plot(x1,f1,'-*','color','k','LineWidth',1.5)
    legend('The FD Line','The Best Fitted Line','Location','NorthEast')
    xlabel('log(r)')
    ylabel('log(N)')
    
    % Figure 3
    f2 = polyval(D2,x2);
    figure, plot(x2,y2,'-*','color','b','LineWidth',1.5)
    grid
    hold on, plot(x2,f2,'-*','color','k','LineWidth',1.5)
    legend('The FD Line','The Best Fitted Line','Location','NorthEast')
    xlabel('log(r)')
    ylabel('log(N)')
end

%% return feature vector
Out = [FD_BC FD_MBC FD_TPSA FDMat_BC FDMat_MBC FDMat_TPSA];

end
