# TODO do we use it?

def findmaxnum(seq):
"""
        this function is obtain the maximum numbers of the given sequence
        note the sequence is stored in cell mode
"""
#
# See also zigzag
#
# Author:
# ---------------------------------------------
#    (C)Xunkai Wei <xunkai.wei@gmail.com>
#    Beijing Aeronautical Technology Research Center
#    Beijing #9203-12,10076
# History:
# ---------------------------------------------
# Creation: beta  Date: 01/10/2007
# Revision: 1.0   Date: 12/11/2007
#
    if type(seq) == list:

        numseq = len(seq)
        maxnum = 1
        for i in range(numseq):
            temp = seq[i]
            numseq = len(temp)
            if numseq > maxnum:
                maxnum = numseq
    else:
        raise ValueError('I was only designed to handle cell sequence')
