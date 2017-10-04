import numpy as np
from ImFEATbox.__helperCommands import conv2float, _complex_dtype
from scipy.misc import factorial

def ZernikeF(I, returnShape=False):
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


    if returnShape:
        return (92,1)

    # reserve space for the variables
    Z = np.zeros(20, dtype=_complex_dtype())
    A = np.zeros(20)
    Phi = np.zeros(20)

    ## determine the moments for different orders
    # TODO ???????????? order always = number of repetition ?
    for m in range(2,42,2): #m=2:2:40
        n = m

        # convert image to square size image
        s = np.shape(I)
        if s[0] < s[1]:
            p = conv2float(I[:,0:s[0]])
        elif s[0] > s[1]:
            p = conv2float(I[0 : s[1],:])
        elif s[0] == s[1]:
            p = conv2float(I)

        # determine Zernike moments
        N = np.shape(p)[0]
        x = range(1, N+1) #1:N
        y = x
        [X,Y] = np.meshgrid(x,y)
        R = np.sqrt(np.power((2*X)-N-1,2)+(np.power((2*Y)-N-1,2)))/float(N)
        Theta = np.arctan2((N-1-(2*Y)+2),((2*X)-N+1-2))
        #R[R<=1] = (R[R<=1] * R[R<=1])

        #R[np.less_equal(R,1)] = R[np.less_equal(R,1)] * R[np.less_equal(R,1)]
        #R = np.reshape(np.less_equal(R,np.ones(R.shape))*1, R.shape ) * R
        R[R>1] = 0
        #R[np.less(R,1)] = 0
        #plt.imshow(R);plt.show()

        # compute Zernike Polynomials
        Rad = np.zeros(R.shape, dtype=_complex_dtype())

        # TODO: for loop does not make sense...?
        for s in range(0,(n-abs(m))/2+1): # 0:(n-abs(m))/2
            c = np.power(-1,s)*factorial(n-s)/(factorial(s)*factorial((n+np.abs(m))/2.0-s)*factorial((n-np.abs(m))/2.0-s))
            Rad = Rad + c * np.power(R,(n-2*s))

        # calculate moments
        Product = p*Rad*np.exp(-1j*m*Theta)
        Zernike = np.sum(Product)

        # count number of pixels inside the unit circle and normalize moments
        cnt = np.count_nonzero(R)+1
        Z[m/2-1] = (n+1)*Zernike / float(cnt)

        # calculate amplitude and phase (in degrees) of the moments
        A[m/2-1] = np.real(np.abs(Zernike))
        Phi[m/2-1] = np.real(np.angle(Zernike)*180/np.pi)

    # calculate mean, std and max values for Z, A and Phi
    mean_Z = np.mean(Z)
    mean_A = np.mean(A)
    mean_Phi = np.mean(Phi)

    std_Z = np.std(Z, ddof=1)
    std_A = np.std(A, ddof=1)
    std_Phi = np.std(Phi, ddof=1)

    max_Z = Z[np.argmax(np.abs(Z))]
    max_A = np.max(A)
    max_Phi = np.max(Phi)

    ## return feature vector

    Out = np.hstack([np.real(Z), np.imag(Z), A, Phi, np.real(mean_Z), np.imag(mean_Z), mean_A, mean_Phi,
        np.real(std_Z), np.imag(std_Z), std_A, std_Phi, np.real(max_Z), np.imag(max_Z), max_A, max_Phi])

    return Out
