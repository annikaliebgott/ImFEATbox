from numpy import genfromtxt
import numpy as np

def DHankelF(I, typeflag):
"""
     Input:     - I: A 2D image
                - typeflag: Struct of logicals to permit extracting features
                  based on desired characteristics:
                       + typeflag.global: all features
                       + typeflag.transform: all features
                       + typeflag.corr: only features based on correlation
                  default: all features are being extracted
                  For more information see README.txt


     Output:    - Out: A (1x75) vector containing 75 metrics calculated from
                       the discrete Hankel transform
"""
    # This implementation uses an array of Bessel roots, which is stored in
    # 'dht.mat'. By default, the first 4097 roots of the Bessel functions of
    # the first kind with order 0 to 4 were precomputed. This allows DHT up to
    # an order of 4 with up to 4096 sampling points by default. Use JNROOTS
    # for calculating more roots.
    #
    # ************************************************************************
    # Modified for MRI feature extraction by the Department of Diagnostic
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
    #
    # Implementation based on:  M. Guizar-Sicairos, J.C. Gutierrez-Vega,
    #                           Computation of quasi-discrete Hankel transforms
    #                           of integer order for propagating optical wave
    #                           fields, J. Opt. Soc. Am. A 21, 53-58 (2004).
    #
    #  Implemented by:  Marcel Leutenegger, June 2006
    #                   Manuel Guizar-Sicairos, 2004


    if 'typeflag' not in globals():
       typeflag.global = True
       typeflag.transform = True
       typeflag.corr = True

    # converte image to float
    I = np.array(I, dtype='float')

    # load precalculated Bessel Jn roots
    c = genfromtxt('bessel.csv', delimiter=',')

    # sizes original image
    N1, N2 = np.shape(I)

    # crop image to a even sized image
    if (N1 % 2) == 0 and (N2 % 2) == 0:
        IMG = I
    if(mod(N1,2)>0)
        IMG = I(1:N1-1,:)
    end
    if(mod(N2,2)>0)
        IMG = I(:,1:N2-1)
    end

    # sizes cropped image
    N4 = size(IMG,2)
    N3 = size(IMG,1)

    # transform order n, signal factor R
    n = 0
    R = 5

    C1=c(n,N1)
    C2=c(n,N2)
    c1=c(n,0:N1)
    c2=c(n,0:N2)
    I1=abs(besselj(1+n,c1))
    I2=abs(besselj(1+n,c2))
    K1=2*pi*R/C1*I1(:)
    K2=2*pi*R/C2*I2(:)
    R1=I1(:)/R
    R2=I2(:)/R
    I1=sqrt(2/C1)./I1
    I2=sqrt(2/C2)./I2
    I1=I1(:)*I1.*besselj(n,c1(:)/C1*c1)
    I2=I2(:)*I2.*besselj(n,c2(:)/C2*c2)
    I2 = I2(1:N4,1:N4)
    I1 = I1(1:N3,1:N3)
    R1 = R1(1:N3)
    R2 = R2(1:N4)
    K1 = K1(1:N3)
    K2 = K2(1:N4)


    ## transform image

    # preallocate hankel transform matrix HTM
    HTM = zeros((size(IMG)))

    # along the vertical dimension
    for j=1 : N3
        HTM(j,:) = (IMG(j,:)*I2*(R2).^(-1))*K2
    end

    # along the horizontal dimension
    for i=1 : N4
        HTM(:,i) = I1*(HTM(:,i).'*((R1).^(-1)))*K1
    end


    ## feature extraction
    # square matrix

    if typeflag.global || typeflag.transform
        a = moment(HTM,2)
        b = moment(HTM,4)
        Out = [mean2(HTM) std2(HTM) std(std(HTM)) max(HTM(:)) min(HTM(:)) mean(a) std(a) mean(b) std(b)]
    end


    ## breake image into four equal sized matrices

    # get the center of the image
    half1 = floor(N3/2)
    half2 = floor(N4/2)
    IMG1 = IMG(1:half1, 1:half2)
    IMG2 = IMG(1:half1,half2+1:N4)
    IMG3 = IMG(half1+1:N3, 1:half2)
    IMG4 = IMG(half1+1:N3, half2+1:N4)

    n = 0
    R = 5

    C3=c(1+n,1+half1)
    C4=c(1+n,1+half2)
    c3=c(1+n,1:half1)
    c4=c(1+n,1:half2)
    I3=abs(besselj(1+n,c3))
    I4=abs(besselj(1+n,c4))
    K3=2*pi*R/C3*I3(:)
    K4=2*pi*R/C4*I4(:)
    R3=I3(:)/R
    R4=I4(:)/R
    I3=sqrt(2/C3)./I3
    I4=sqrt(2/C4)./I4
    I3=I3(:)*I3.*besselj(n,c3(:)/C3*c3)
    I4=I4(:)*I4.*besselj(n,c4(:)/C4*c4)


    ## transform decomposed image
    HTM1 = zeros(half1, half2)
    HTM2 = zeros(half1, half2)
    HTM3 = zeros(half1, half2)
    HTM4 = zeros(half1, half2)

    # along the vertical dimension
    for j=1 : half1
        HTM1(j,1:half2) = I4*(IMG1(j,1:half2)*(R4).^(-1))*K4
        HTM2(j,1:half2) = I4*(IMG2(j,1:half2)*(R4).^(-1))*K4
        HTM3(j,1:half2) = I4*(IMG3(j,1:half2)*(R4).^(-1))*K4
        HTM4(j,1:half2) = I4*(IMG4(j,1:half2)*(R4).^(-1))*K4
    end

    # along the horizontal dimension
    for i=1 : size(IMG1,2)
        HTM1(1:half1,i) = I3*(HTM1(1:half1,i).'*((R3).^(-1)))*K3
        HTM2(1:half1,i) = I3*(HTM2(1:half1,i).'*((R3).^(-1)))*K3
        HTM3(1:half1,i) = I3*(HTM3(1:half1,i).'*((R3).^(-1)))*K3
        HTM4(1:half1,i) = I3*(HTM4(1:half1,i).'*((R3).^(-1)))*K3
    end

    ## feature extraction
    # 2nd and 4th moments
    if typeflag.global || typeflag.transform
        a1 = moment(HTM1,2)
        b1 = moment(HTM1,4)
        a2 = moment(HTM2,2)
        b2 = moment(HTM2,4)
        a3 = moment(HTM3,2)
        b3 = moment(HTM3,4)
        a4 = moment(HTM4,2)
        b4 = moment(HTM4,4)

        Out = [Out mean2(HTM1) std2(HTM1) std(std(HTM1)) mean(a1) std(a1) mean(b1) std(b1) max(HTM1(:)) min(HTM1(:)) mean2(HTM2) std2(HTM2) std(std(HTM2)) mean(a2) std(a2) mean(b2) std(b2) max(HTM2(:)) min(HTM2(:)) mean2(HTM3) std2(HTM3) std(std(HTM3)) mean(a3) std(a3) mean(b3) std(b3) max(HTM3(:)) min(HTM3(:)) mean2(HTM4) std2(HTM4) std(std(HTM4)) mean(a4) std(a4) mean(b4) std(b4) max(HTM4(:)) min(HTM4(:))]
    end

    # correlation
    corr12 = corr(HTM1(1:half1,1:half2),HTM2(1:half1,1:half2))
    corr13 = corr(HTM1(1:half1,1:half2),HTM3(1:half1,1:half2))
    corr14 = corr(HTM1(1:half1,1:half2),HTM4(1:half1,1:half2))
    corr32 = corr(HTM3(1:half1,1:half2),HTM2(1:half1,1:half2))
    corr42 = corr(HTM4(1:half1,1:half2),HTM2(1:half1,1:half2))
    corr34 = corr(HTM4(1:half1,1:half2),HTM3(1:half1,1:half2))


    ## return feature vector

    if ~(typeflag.transform || typeflag.global)
        # extract only correlation based features
        Out = [std2(corr12) mean2(corr12) max(corr12(:)) min(corr12(:)) nnz(corr12) std2(corr13) mean2(corr13) max(corr13(:)) min(corr13(:)) nnz(corr13) std2(corr14) mean2(corr14) max(corr14(:)) min(corr14(:)) nnz(corr14) std2(corr32) mean2(corr32) max(corr32(:)) min(corr32(:)) nnz(corr32) std2(corr42) mean2(corr42) max(corr42(:)) min(corr42(:)) nnz(corr42) std2(corr34) mean2(corr34) max(corr34(:)) min(corr34(:)) nnz(corr34)]
    else
        # extract all features
        Out = [Out std2(corr12) mean2(corr12) max(corr12(:)) min(corr12(:)) nnz(corr12) std2(corr13) mean2(corr13) max(corr13(:)) min(corr13(:)) nnz(corr13) std2(corr14) mean2(corr14) max(corr14(:)) min(corr14(:)) nnz(corr14) std2(corr32) mean2(corr32) max(corr32(:)) min(corr32(:)) nnz(corr32) std2(corr42) mean2(corr42) max(corr42(:)) min(corr42(:)) nnz(corr42) std2(corr34) mean2(corr34) max(corr34(:)) min(corr34(:)) nnz(corr34)]
    end

    return Out
