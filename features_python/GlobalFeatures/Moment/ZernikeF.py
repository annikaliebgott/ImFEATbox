import numpy as np

def ZernikeF(I, test=False):
"""
     Input:     - I: A 2D image


     Output:    - Out: A (1x92) vector containing 92 metrics based on Zernike
                  moments
"""
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
    # Copyright C 2014 Amir Tahmasbi
    # Texas A&M University
    # amir.tahmasbi@tamu.edu
    # http://people.tamu.edu/~amir.tahmasbi/index.html
    #
    # License Agreement: To acknowledge the use of the code please cite the
    #                    following papers:
    #
    # [1] A. Tahmasbi, F. Saki, S. B. Shokouhi,
    #     Classification of Benign and Malignant Masses Based on Zernike Moments,
    #     Comput. Biol. Med., vol. 41, no. 8, pp. 726-735, 2011.
    #
    # [2] F. Saki, A. Tahmasbi, H. Soltanian-Zadeh, S. B. Shokouhi,
    #     Fast opposite weight learning rules with application in breast cancer
    #     diagnosis, Comput. Biol. Med., vol. 43, no. 1, pp. 32-41, 2013.



    # reserve space for the variables
    Z = np.zeros(1,20)
    A = np.zeros(1,20)
    Phi = np.zeros(1,20)

    ## determine the moments for different orders
    # TODO ???????????? wat
    for m=2:2:40
        n = m

        # converte image to square size image
        s = size(I)
        if (s(1,1) < s(1,2))
            p = double( I(:,1:s(1,1)) )
        elseif (s(1,1) > s(1,2))
            p = double( I(1 : s(1,2), :) )
        elseif (s(1,1) == s(1,2))
            p = double( I )
        end

        # determine Zernike moments
        N = size(p,1)
        x = 1:N
        y = x
        [X,Y] = meshgrid(x,y)
        R = double( sqrt((2.*X-N-1).^2+(2.*Y-N-1).^2)/N )
        Theta = atan2((N-1-2.*Y+2),(2.*X-N+1-2))
        R = (R<=1).*R

        # compute Zernike Polynomials
        Rad = zeros(size(R))
        for s = 0:(n-abs(m))/2
            c = (-1)^s*factorial(n-s)/(factorial(s)*factorial((n+abs(m))/2-s)*...
                factorial((n-abs(m))/2-s))
            Rad = Rad + c*R.^(n-2*s)
        end
        Rad = double(Rad)

        # calculate moments
        Product = p(x,y).*Rad.*exp(-1i*m*Theta)
        Zernike = sum(Product(:))

        # count number of pixels inside the unit circle and normalize moments
        cnt = nnz(R)+1
        Z(m/2) = (n+1)*Zernike/cnt

        # calculate amplitude and phase (in degrees) of the moments
        A(m/2) = real(double(abs(Zernike)))
        Phi(m/2) = real(double(angle(Zernike)*180/pi))

        # calculate mean, std and max values for Z, A and Phi
        mean_Z = mean(Z)
        mean_A = mean(A)
        mean_Phi = mean(Phi)

        std_Z = std(Z)
        std_A = std(A)
        std_Phi = std(Phi)

        max_Z = max(Z)
        max_A = max(A)
        max_Phi = max(Phi)
    end


    ## return feature vector

    Out = [real(Z) imag(Z) A Phi real(mean_Z) imag(mean_Z) mean_A mean_Phi...
        real(std_Z) imag(std_Z) std_A std_Phi real(max_Z) imag(max_Z) max_A max_Phi]

    return Out
