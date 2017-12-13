import numpy as np

def rle_0(si, NL):
    """
     RLE   image gray level Run Length matrix for 0 degree
    """
    #
    # Author:
    # ---------------------------------------------
    #    (C)Xunkai Wei <xunkai.wei@gmail.com>
    #    Beijing Aeronautical Technology Research Center
    #    Beijing #9203-12,10076
    # History:
    #  -------
    # Creation: beta  Date: 01/11/2007
    # Revision: 1.0   Date: 10/11/2007

    # Assure row number is exactly the gray level
    m, n = np.shape(si)

    oneglrlm = np.zeros((NL, n))

    print("NL=" + str(NL) + ", n=" + str(n))

    for i in range(m):
        x = si[i,:]
        # run length Encode of each vector
        index = np.where(x[:-1] != x[1:])[0][:len(x)]
        # run lengths
        lenX = np.diff(np.append(0, index+1))
        # run values
        print(x.min())
        val = x[index]
        # compute current numbers (or contribution) for each bin in GLRLM
        temp = np.zeros((NL,n))
        temp[val,lenX] = 1
        #temp = np.ufunc.at(np.hstack([val, lenX]), 1, np.hstack([NL, n]))
        oneglrlm = temp + oneglrlm # accumulate each contribution
    return oneglrlm
