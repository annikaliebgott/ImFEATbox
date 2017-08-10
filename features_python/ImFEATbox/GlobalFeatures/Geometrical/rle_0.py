def rle_0(si, NL):
"""
     RLE   image gray level Run Length matrix for 0degree
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

    oneglrlm = zeros(NL, n)

    for i in range(m):
        x = si[i,:]
        # run length Encode of each vector
        index = [ np.where(x[:-2] != x[1:]), len(x) ]
        lenX = np.diff([ 0, index ]) # run lengths
        val = x[index]          # run values
        temp = accumarray([val, lenX].T, 1, [NLn n]) # compute current numbers (or contribution) for each bin in GLRLM
        oneglrlm = temp + oneglrlm # accumulate each contribution
    return oneglrlm
