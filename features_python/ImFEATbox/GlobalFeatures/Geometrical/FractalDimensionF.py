from math import log
from PIL import Image
from scipy.misc import imresize
import numpy as np

def FractalDimensionF(I, width=256, plotflag=False, test=False):
"""
 Input:     - I: A 2D image
            - plotflag: A locical flag to enable/disable visualization.
              Default: False
            - width: largest size of the box. Default: width = 256
            - test: A logical flag to return only the shape of features.
              Default: False

 Output:    - Out: A (1 x (log(width)/log(2))*3 + 3) vector containing
                   metrics based on self-similarity of image structures


 Implemented algorithms:   1. Box-Counting(Haus) (BC)
                           2. Differential Box-Counting (MBC)
                           3. Triangular Prism Surface Area (TPSA)
"""
# ************************************************************************
# Implemented for MRI feature extraction by the Department of Diagnostic
# and Interventional Radiology, University Hospital of Tuebingen, Germany
# and the Institute of Signal Processing and System Theory University of
# Stuttgart, Germany. Last modified: November 2016
#
# This implementation is part of ImFEATbox, a toolbox for image feature
# extraction and analysis. Available online at:
# https://github.com/annikaliebgott/ImFEATbox
#
# Contact: annika.liebgott@iss.uni-stuttgart.de
# ************************************************************************

    #if ~exist('plotflag','var')
    #    plotflag = false;
    #end
    #if ~exist('width','var')
    #    width = 256;
    #end

    if I.ndim > 2:
        raise ValueError('Maximum dimension should be 2.')

    # Spliting the Image into the small grids and copute the fractal dimension

    # The largest size of the box
    p = int(log(width)/log(2))
    RescaledI = imresize(I, (width, width), interp='bicubic')

    # Allocation of the number of box size
    counter = 0
    counter_dbc = 0
    counter_tpsa = 0
    step = width/np.power(2,range(1,p+1))
    N_mbc = np.zeros(len(step))
    N_tpsa = np.zeros(len(step))
    N_b = np.zeros(len(step))

    ## 2D boxcount
    for n in range(0,len(step)):
        stepnum = step(n);
        for i in range(0, width, stepnum):
            for j in range(0, width, stepnum):

                # Get the Grid in each level
                testim = RescaledI[i:i+stepnum-1,j:j+stepnum-1]

                # Differential(Modified) Box Counting
                MaxGrayLevel = np.max(testim[:])
                MinGrayLevel = np.min(testim[:])
                GridCont = MaxGrayLevel-MinGrayLevel+1
                counter_dbc = counter_dbc + GridCont
                # Differential(Modified) Box Counting (MBC)

                #Triangular Prism Surface Area (TPSA)
                a = testim[1,1]
                b = testim[1,-1]
                c = testim[-1,1]
                d = testim[-1,-1]
                e = (a+b+c+d)/4

                w = sqrt(np.power((b-a),2) + np.power(stepnum,2))
                x = sqrt(np.power((c-b),2) + np.power(stepnum,2))
                y = sqrt(np.power((d-c),2) + np.power(stepnum,2))
                z = sqrt(np.power((a-d),2) + np.power(stepnum,2))

                o = sqrt(np.power((a-e),2) + np.power(0.5*stepnum,2))
                p2 = sqrt(np.power((b-e),2) + np.power(0.5*stepnum,2))
                q = sqrt(np.power((c-e),2) + np.power(0.5*stepnum,2))
                t = sqrt(np.power((d-e),2) + np.power(0.5*stepnum,2))

                # Using Herons Formula
                sa = (w+p2+o)/2
                sb = (x+p2+q)/2
                sc = (y+q+t)/2
                sd = (z+o+t)/2

                # Areas of Traiangle
                S_ABE = sqrt(sa*(sa-w)*(sa-p2)*(sa-o))
                S_BCE = sqrt(sb*(sb-x)*(sb-p2)*(sb-q))
                S_CDE = sqrt(sc*(sc-q)*(sc-t)*(sc-y))
                S_DAE = sqrt(sd*(sd-z)*(sd-o)*(sd-t))
                SurfaceArea = S_ABE + S_BCE + S_CDE + S_DAE
                counter_tpsa = counter_tpsa + SurfaceArea
                #Triangular Prism Surface Area


                # Basic Box Counting (BC)
                #if (size(find(testim~=0),1)~=0)
                # TODO test this
                if np.size(np.nonzero(testim)) > 0:
                    counter += 1

        N_mbc[n] = counter_dbc
        N_tpsa[n] = counter_tpsa
        N_b[n] = counter
        counter = 0
        counter_dbc = 0
        counter_tpsa = 0


    # Resolution
    r0 = np.power(2,(range(p,0,-1)))

    # Dimension of BC
    x0 = np.log(r0)
    y0 = np.log(N_b)
    FDMat_BC = y0/x0
    D0 = np.polyfit(x0, y0, 1)
    FD_BC = D0[0]

    # Dimension of MBC
    x1 = np.log(r0)
    y1 = np.log(N_mbc)
    FDMat_MBC = y1/x1
    D1 = np.polyfit(x1, y1, 1)
    FD_MBC = D1[0]

    # Dimension of TPSA
    x2 = log(r0)
    y2 = log(N_tpsa)
    FDMat_TPSA = y2/x2
    D2 = np.polyfit(x2, y2, 1)
    FD_TPSA = 2 - D2[0]

    # Plotting
    if plotflag:
        import matplotlib.pyplot as plt

        # Figure 1
        f0 = np.polyval(D0,x0)
        plt.plot(x0,y0,'-*', color='b', linewidth=1.5, label='The FD Line')
        plt.grid()
        plt.plot(x0,f0,'-*','color','k', linewidth=1.5, label='The Best Fitted Line')
        legend()
        plt.xlabel('log(r)')
        plt.ylabel('log(N)')

        # Figure 2
        plt.figure()
        f1 = np.polyval(D1,x1)
        plt.plot(x1,y1,'-*', color='b', linewidth=1.5, label='The FD Line')
        plt.grid()
        plt.plot(x1,f1,'-*','color','k', linewidth=1.5, label='The Best Fitted Line')
        legend()
        plt.xlabel('log(r)')
        plt.ylabel('log(N)')

        # Figure 3
        plt.figure()
        f2 = np.polyval(D2,x2)
        plt.plot(x2,y2,'-*','color','b', linewidth=1.5, label='The FD Line')
        plt.grid()
        plt.plot(x2,f2,'-*','color','k', linewidth=1.5, label='The Best Fitted Line')
        legend()
        plt.xlabel('log(r)')
        plt.ylabel('log(N)')

        plt.show()


    ## return feature vector
    return [FD_BC, FD_MBC, FD_TPSA, FDMat_BC, FDMat_MBC, FDMat_TPSA]
