function Out = QuadtreeDecompositionF(I,threshold)
% Input:     - I: A 2D image
%            - threshold: threshold for the quadtree decomposition
%              (default: threshold = 0.27)
% 
%
% Output:    - Out: A (1x5) vector containing 5 metrics calculated from
%              the quadtree decomposition of the image
%
% ************************************************************************
% Implemented for MRI feature extraction by the Department of Diagnostic 
% and Interventional Radiology, University Hospital of Tuebingen, Germany 
% and the Institute of Signal Processing and System Theory University of 
% Stuttgart, Germany. Last modified: November 2016
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
% ************************************************************************

if ~exist('threshold', 'var')
    threshold = 0.27;
end    

B = double(I);

% extend I to a square matrix of the next power of 2 by zero padding
s = size(B);    
if (s(1,1) < s(1,2))
    B = padarray(B, [2^( nextpow2( s(1,2) )) - s(1,1) ; 2^(nextpow2(s(1,2)))-s(1,2)],'post');
elseif (s(1,1) > s(1,2))
    B = padarray(B, [2^(nextpow2(s(1,1)))- s(1,1); 2^(nextpow2(s(1,1)))- s(1,2)],'post');
elseif (s(1,1) == s(1,2))
    B = padarray(B, [ 2^(nextpow2(s(1,1)))- s(1,1); 2^(nextpow2(s(1,1)))- s(1,2)], 'post');
end


%% decompose image

S = qtdecomp(B,threshold);    
blocks = repmat(uint8(0),size(S));
sum_blocks = 0;
mean_r = 0;
mean_c = 0;

%% extract features
for dim = [2048 1024 512 256 128 64 32 16 8 4 2 1];
    numblocks = length(find(S==dim));
    if (numblocks > 0)
        values = repmat(uint8(1),[dim dim numblocks]);
        values(2:dim,2:dim,:) = 0;
        blocks = qtsetblk(blocks,S,dim,values);
        sum_blocks = 1 + sum_blocks;
        
        [~,r,c] = qtgetblk(blocks,S,dim);   
        mean_r = mean2(r)+ mean_r;
        mean_c = mean2(c)+ mean_c;
        
    end
end

blocks(end,1:end) = 1;
blocks(1:end,end) = 1;

% length
length_quadtree = nnz(blocks);

% number of non zero elements
entries = nnz(S);


%% return feature vector
Out = [sum_blocks length_quadtree mean_r mean_c entries];

end    