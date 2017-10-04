## Function to compute the Run length matrix

# The script has been written by Xunkai Wei <xunkai.wei@gmail.com>
# Beijing Aeronautical Technology Research Center

def grayrlmatrix(I, Offset=np.array([1,2,3,4]), NumLevels=None, GrayLimits=None):
"""
  Description
  -------------------------------------------
   Computes the graylevel run length (GLRL) matrix used for textural
   analysis of an image using zigzag scan method.The method includes four
   basic steps
       Step 1 determine direction
       Step 2 zigzag scan
       Step 3 obtain new sequences
       Step 4 calculate run-length matrix
   -----------------------------------------
   GLRLMS = GRAYRLMATRIX(I,PARAM1,VALUE1,PARAM2,VALUE2,...) returns one or more
   gray-level run-length matrices, depending on the values of the optional
   parameter/value pairs. Parameter names can be abbreviated, and case does
   not matter.
  ------------------------------------------
   Parameters include:
  ------------------------------------------
   'Offset'         A p-by-1 vector of offsets specifying the scanning direction.


                    Angle     OFFSET
                    -----     ------
                    0          1
                    45         2
                    90         3
                    135        4

                    OFFSET must be integers from {1 2 3 4}.

                    Default: [1 2 3 4]

   'NumLevels'      An integer specifying the number of gray levels to use when
                    scaling the grayscale values in I. For example, if
                    'NumLevels' is 8, GRAYRLMATRIX scales the values in I so
                    they are integers between 1 and 8.  The number of gray levels
                    determines the size of the gray-level run-length matrix


                    'NumLevels' must be an integer. 'NumLevels' must be 2 if I
                    is logical.

                    Default: 8 for numeric
                             2 for logical

   'GrayLimits'     A two-element vector, [LOW HIGH], that specifies how the
                    grayscale values in I are linearly scaled into gray
                    levels. Grayscale values less than or equal to LOW are
                    scaled to 1. Grayscale values greater than or equal to
                    HIGH are scaled to HIGH.  If 'GrayLimits' is set to [],
                    GRAYRLMATRIX uses the minimum and maximum grayscale values
                    in I as limits, [min(I(:)) max(I(:))].

                    Default: the LOW and HIGH values specified by the
                    class, e.g., [LOW HIGH] is [0 1] if I is double and
                    [-32768 32767] if I is int16.

  ------------------------------------------
  Example
  ------------------------------------------
 I =[1     1     1     2     2
      3     4     2     2     3
      4     4     4     4     4
      5     5     3     3     3
      1     1     3     4     5]
 [GLRLMS,SI] = grayrlmatrix(I,'NumLevels',5,'G',[])
 I =
      1     1     1     2     2
      3     4     2     2     3
      4     4     4     4     4
      5     5     3     3     3
      1     1     3     4     5
 GLRLMS(:,:,1) =
      0     1     1     0     0
      0     2     0     0     0
      3     0     1     0     0
      2     0     0     0     1
      1     1     0     0     0
 GLRLMS(:,:,2) =
      5     0     0     0     0
      0     2     0     0     0
      4     1     0     0     0
      5     1     0     0     0
      3     0     0     0     0
 GLRLMS(:,:,3) =
      5     0     0     0     0
      2     1     0     0     0
      4     1     0     0     0
      5     1     0     0     0
      3     0     0     0     0
 GLRLMS(:,:,4) =
      5     0     0     0     0
      4     0     0     0     0
      6     0     0     0     0
      5     1     0     0     0
      3     0     0     0     0
 SI =
      1     1     1     2     2
      3     4     2     2     3
      4     4     4     4     4
      5     5     3     3     3
      1     1     3     4     5
 -------------------------------------------
 See also zigzag rle_0 rle_45
 -------------------------------------------
"""
# Author:
# -------------------------------------------
#    (C)Xunkai Wei <xunkai.wei@gmail.com>
#    Beijing Aeronautical Technology Research Center
#    Beijing #9203-12,100076
# -------------------------------------------
# History:
# -------------------------------------------
# Creation: beta  Date: 01/10/2007
# Revision: 1.0   Date: 14/11/2007
# -------------------------------------------
# Bug Fixed:
# -------------------------------------------
# 1.Issue wrong results for nonsquare matrix,now output cells instead of
#   multi-dim arrays
# 2.Add support for inputs checking inspired by MATLAB style
#
    [I, Offset, NumLevels, GrayLimits] = ParseInputs(varargin{:})

    if NumLevels == None:
        if I.dtype == bool:
            NumLevels = 2
        else:
            NumLevels = 8

    # gl = getrangefromclass(I)
    if GrayLimits == None:
        if I.dtype == bool:
            GrayLimits = [0 , 1]
        elif I.dtype == 'uint8':
            GrayLimits = [0 , 255]
        elif I.dtype == float:
            GrayLimits = [0 , 1]



    # Scale I so that it contains integers between 1 and NumLevels.
    if GrayLimits[1] == GrayLimits[0]:
        SI = np.ones(np.shape(I))
    else:
        slope = (NumLevels-1) / (GrayLimits[1] - GrayLimits[0])
        intercept = 1 - (slope*(GrayLimits[0])
        SI = np.round(imlincomb(slope,I,intercept,'double'))

    # Clip values if user had a value that is outside of the range, e.g., double
    # image = [0 .5 2;0 1 1]; 2 is outside of [0,1]. The order of the following
    # lines matters in the event that NumLevels = 0.
    SI[SI > NumLevels] = NumLevels
    SI[SI < 1] = 1
    # total numbers of direction
    numOffsets = np.size(Offset,1)

    GLRLMS = np.array([])
    if NumLevels != 0:
        # make direction matrix for all given directions
        for k in range(numOffsets):
            np.append(GLRLMS, computeGLRLM(SI,Offset[k],NumLevels))


# --------------------------------------------------------------------

def computeGLRLM(si,offset,NumLevels):
    # For given direction, compute the run length matrix

    if offset == 1:
        # 0 degree
        oneGLRLM = rle_0(si,NumLevels)
    if offset == 2:
        # 45 degree
        seq = zigzag(si)
        oneGLRLM  = rle_45(seq,NumLevels)
    if offset == 3:
        # 90 degree
        oneGLRLM = rle_0(si.T,NumLevels)
    if offset == 4:
        # 135 degree
        seq = zigzag(fliplr(si))
        oneGLRLM = rle_45(seq,NumLevels)
    else:
        raise ValueError('Only 4 directions supported')

    return oneGLRLM


# TODO end here
    # --------------------------------------------------------------------
    function [I, offset, NumLevels, gl] = ParseInputs(varargin)
    # parsing parameter checking
    # Inputs must be max seven item
    iptchecknargin(1,7,nargin,mfilename)
    #
    # Check I
    I = varargin{1}
    iptcheckinput(I,{'logical','numeric'},{'2d','real','nonsparse'}, ...
        mfilename,'I',1)
    # ------------------------
    # Assign Defaults
    # -------------------------
    # four directions 0, 45, 90,135
    offset = [1;2;3;4]

    gl = getrangefromclass(I)

    # Parse Input Arguments
    if nargin ~= 1

        paramStrings = {'Offset','NumLevels','GrayLimits'}

        for k = 2:2:nargin

            param = lower(varargin{k})
            inputStr = iptcheckstrs(param, paramStrings, mfilename, 'PARAM', k)
            idx = k + 1  #Advance index to the VALUE portion of the input.
            if idx > nargin
                eid = sprintf('Images:#s:missingParameterValue', mfilename)
                msg = sprintf('Parameter ''#s'' must be followed by a value.', inputStr)
                error(eid,'#s', msg)
            end

            switch (inputStr)

                case 'Offset'

                    offset = varargin{idx}
                    iptcheckinput(offset,{'logical','numeric'},...
                        {'d','nonempty','integer','real'},...
                        mfilename, 'OFFSET', idx)
                    # must be row vector
                    if size(offset,2) ~= 1
                        eid = sprintf('Images:#s:invalidOffsetSize',mfilename)
                        msg = 'OFFSET must be an n x 1 array.'
                        error(eid,'#s',msg)
                    end
                    offset = double(offset)

                case 'NumLevels'

                    NumLevels = varargin{idx}
                    iptcheckinput(NumLevels,{'logical','numeric'},...
                        {'real','integer','nonnegative','nonempty','nonsparse'},...
                        mfilename, 'NumLevels', idx)
                    if np.prod(np.shape(NumLevels)) > 1:
                        eid = sprintf('Images:#s:invalidNumLevels',mfilename)
                        msg = 'NumLevels cannot contain more than one element.'
                        error(eid,'#s',msg)
                    elseif islogical(I) && NumLevels ~= 2
                        eid = sprintf('Images:#s:invalidNumLevelsForBinary',mfilename)
                        msg = 'NumLevels must be two for a binary image.'
                        error(eid,'#s',msg)
                    end
                    NumLevels = double(NumLevels)

                case 'GrayLimits'

                    gl = varargin{idx}
                    iptcheckinput(gl,{'logical','numeric'},{'vector','real'},...
                        mfilename, 'GL', idx)
                    if isempty(gl)
                        gl = [min(I(:)) max(I(:))]
                    elif np.prod(np.shape(gl)) != 2:
                        eid = sprintf('Images:#s:invalidGrayLimitsSize',mfilename)
                        msg = 'GL must be a two-element vector.'
                        error(eid,'#s',msg)
                    end
                    gl = double(gl)
            end
        end
    return GLRLMS, SI
