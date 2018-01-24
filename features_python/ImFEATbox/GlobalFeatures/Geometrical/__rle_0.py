import numpy as np
#from ImFEATbox.__helperCommands import *
from ImFEATbox.__helperCommands import _int_dtype

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

    oneglrlm = np.zeros((NL-1, n))

    #print("NL=" + str(NL) + ", n=" + str(n))

    for i in range(m):
        x = si[i,:]
        # run length Encode of each vector
        index = np.append(np.where(x[:-1] != x[1:])[0], len(x)-1)
        # run lengths
        lenX = np.array(np.diff(np.append(0, index+1)), dtype=_int_dtype())
        # run values
        #print(x.min())
        val = np.array(x[index], dtype=_int_dtype())
        # compute current numbers (or contribution) for each bin in GLRLM
        temp = np.zeros((NL-1,n))
        #temp[val,lenX] = 1
        tmp=np.array([val , lenX]).T


        # TODO: count occurrence of (x,y) in tmp and save it in matrix at position (x,y) here
        for j in range(tmp.shape[0]):
            temp[tmp[j]] += 1
        print(temp[:12,:12])

        #temp = np.ufunc.at(np.hstack([val, lenX]), 1, np.hstack([NL, n]))
        oneglrlm = temp + oneglrlm # accumulate each contribution
    return oneglrlm
