import numpy as np
import skimage
from scipy.ndimage.filters import rank_filter

def findLocalMaximum(val, radius):
    """
     Determine the local maximum of a given value

     Author :: Vincent Garcia
     Date   :: 09/02/2007

     INPUT
     =====
     val    : the NxM matrix containing values
     radius : the radius of the neighborhood

     OUTPUT
     ======
     row       : the row position of the local maxima
     col       : the column position of the local maxima
     max_local	: the NxM matrix containing values of val on unique local maximum

     EXAMPLE
     =======
     [l,c,m] = findLocalMaximum(img,radius)
    """

    # FIND UNIQUE LOCAL MAXIMA USING FILTERING (FAST)
    mask = skimage.morphology.disk(r)

    nb    = np.sum(mask)
    highest = rank_filter(val, nb, footprint=mask==1)
    second_highest   = rank_filter(val, nb-1, mask)
    index            = highest==val & highest !=second_highest
    max_local        = zeros(size(val))
    max_local(index) = val(index)
    [row,col]        = find(index==1)


return [row, col, max_local]
