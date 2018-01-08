import numpy as np

def zigzag(SI):
    #
    #  Description:
    #  ------------
    #  This function is used to build the corresponding sequences of a given
    #  scaled gray level image matrix from 45' degree direction. The whole process is using zigzag method
    #  It can handle nonsquare image matrix
    #
    # Author:
    # -------
    #    (C)Xunkai Wei <xunkai.wei@gmail.com>
    #    Beijing Aeronautical Technology Research Center
    #    Beijing #9203-12,10076
    #
    # History:
    #  -------
    # Creation: beta  Date: 01/11/2007
    # Revision: 1.0   Date: 12/11/2007
    #
    # Trick: all the sequence starts or ends lie on the boundary.

    # initializing the variables
    #----------------------------------
    c = 1 # initialize colum indicator
    r = 1 # initialize row   indicator

    rmin = 1 # row   boundary checker
    cmin = 1 # colum boundary checker

    rmax, cmax = np.shape(SI) # get row and colum numbers

    i = 1 # counter for current ith element
    j = 0 # indicator for determining sequence interval

    seq = []

    # intialize sequence mark
    sq_up_begin = 1

    sq_down_begin = 1

    # # Output contain value and its flag status
    # the first row contain value
    # the second row contain its flag
    output = np.zeros( (rmax) * (cmax) )
    # sequence counter
    #
    # # Position Matrix
    # position =zeros(1, rmax * cmax)
    #----------------------------------

    while (r <= rmax) and (c <= cmax):
        # for current point, judge its zigzag direction up 45, or down 45, or
        # 0,or down 90

        if (c + r % 2) == 0:     # up 45 direction
            #  if we currently walk to the left first colum
            if r == rmin:
                # First, record current point
                output[i] = SI[r-1, c-1]
                # if we walk to right last colum
                if c == cmax:
                    # add row number move straight down 90
                    r = r + 1
                    sq_up_end = i
                    sq_down_begin = i+1
                    if len(seq) <= j:
                        seq.append([])
                    seq[j] = output[sq_up_begin-1:sq_up_end-1]
                    j = j + 1
                else:
                    # Continue to move to next (1,c+1) point
                    # This next point should be the begin point of next sequence
                    c = c + 1
                    sq_up_end = i
                    sq_down_begin = i + 1
                    if len(seq) <= j:
                        seq.append([])
                    seq[j] = output[sq_up_begin-1:sq_up_end-1]
                    j = j + 1
                # add couter
                i = i + 1
                # if we currently walk to the last column
            elif (c == cmax) and (r < rmax):
                # first record the point
                output[i] = SI[r-1, c-1]
                # then move straight down to next row
                r = r + 1

                sq_up_end = i
                if len(seq) <= j:
                    seq.append([])
                seq[j] = output[sq_up_begin-1:sq_up_end-1]
                sq_down_begin = i + 1
                j = j + 1
                # add counter
                i = i + 1
                # all other cases i.e. nonboundary points
            elif (r > rmin) and (c < cmax):
                output[i] = SI[r-1, c-1]
                # move to next up 45 point
                r = r - 1
                c = c + 1
                # add counter
                i = i + 1
            # down 45 direction
        else:
            # if we walk to the last row
            if (r == rmax) and (c <= cmax):
                # firstly record current point
                output[i] = SI[r-1, c-1]
                # move right to next point
                c = c + 1
                sq_down_end = i
                if len(seq) <= j:
                    seq.append([])
                seq[j] = output[sq_down_begin-1:sq_down_end-1]
                sq_up_begin = i + 1
                j = j + 1
                # add counter
                i = i + 1
                # if we walk to the first column
            elif c == cmin:
                #first record current point
                output[i] = SI[r-1, c-1]
                if r == rmax:
                    c = c + 1
                    sq_down_end = i
                    if len(seq) <= j:
                        seq.append([])
                    seq[j] = output[sq_down_begin-1:sq_down_end-1]
                    sq_up_begin = i + 1
                    j = j + 1
                else:
                    r = r + 1
                    # record sequence end
                    sq_down_end = i
                    if len(seq) <= j:
                        seq.append([])
                    seq[j] = output[sq_down_begin-1:sq_down_end-1]
                    sq_up_begin =i+1
                    j = j + 1

                i = i + 1
                # all other cases without boundary point
            elif (r < rmax) and (c > cmin):
                output[i] = SI[r-1, c-1]
                # position(i) = sub2ind(SI,r,c)
                r = r + 1
                c = c - 1
                # keep down_info
                i = i + 1

        if (r == rmax) and (c == cmax): # bottom right element
            output[i] = SI[r-1, c-1]
            sq_end = i
            if len(seq) <= j:
                seq.append([])
            seq[j] = output[sq_end-1]
            # position(i) = sub2ind(SI,r,c)

            break
    #seq = np.hstack(seq)
    return seq
