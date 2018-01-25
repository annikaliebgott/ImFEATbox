import numpy as np
from ImFEATbox.__helperCommands import _int_dtype

def rle_45(seq, NL):
    """
         RLE   image gray level Run Length matrix for 45 and 135
         This file is to handle the zigzag scanned sequence for 45 or 135 degree
         direction. Note for 135, just swap the left and the right colum
         seq: zigzag scanned Sequence
         NL: NumLevels
    """
    # Author:
    # ---------------------------------------------
    #    (C)Xunkai Wei <xunkai.wei@gmail.com>
    #    Beijing Aeronautical Technology Research Center
    #    Beijing #9203-12,10076
    # History:
    #  -------
    # Creation: beta  Date: 01/11/2007
    # Revision: 1.0   Date: 10/11/2007


    # Assure row number is exactly the gray leve;
    # number of seqence
    #m = len(seq)

    # number to store the possible max coloums
    #n = np.argmax(np.array(seq))
    n = 1

    for s in range(len(seq)):
        #print(np.shape(seq))
        #print(seq[s])
        if np.size(seq[s]) > 1:
            n = max(n,len(seq[s]))
    #print("n=" + str(n))

    oneglrlm = np.zeros((NL-1, n), dtype=_int_dtype())

    for i in range(len(seq)):
        x = seq[i]

        # run length Encode of each vector
        #index = np.where(x[:-1] != x[1:])[0][:len(x)]
        #index = np.where(x[0:] != x[1:])[0][:len(x)] # TODO: check if same output as matlab
        if np.size(seq[i]) > 1:
            index = np.append(np.where(x[:-1] != x[1:])[0], len(x)-1)
            val = np.array(x[index], dtype=_int_dtype())
        #index = [ find(x(1:end-1) != x(2:end)), length(x) ]
        else:
            index = 0
            val = x
        lenX = np.array(np.diff(np.append(0, index+1)), dtype=_int_dtype())
        # run values
        #print(x.min())

        temp = np.zeros((NL-1,n), dtype=_int_dtype())
        tmp=np.array([val , lenX]).T - 1
        #print("val " + str(val))
        #print("lenX " + str(lenX))
        #if val != [] and lenX != []:

        if tmp.shape == (2,):
            temp[tmp[0],tmp[1]] += 1
            #temp = accumarray([val;len].T, 1, [NL, n]) # compute current numbers (or contribution) for each bin in GLRLM
        else:
            for j in range(tmp.shape[0]):
                temp[tmp[j,0],tmp[j,1]] += 1
                #temp = accumarray([val;len].T, 1, [NL, n]) # compute current numbers (or contribution) for each bin in GLRLM
        oneglrlm = temp + oneglrlm # accumulate each contribution

    return oneglrlm
