function [ SegmentedImage, SegmentedBackground ] = fSegmentCrop ( I, BinaryImage, cropMode )
% segment and crop the image

if(ismatrix(BinaryImage))
    [iRow, iCol] = find(BinaryImage==1);
    iSli = 1;
    iLin = sub2ind(size(BinaryImage), iRow, iCol);
elseif(ndims(BinaryImage) == 3)
    iLin = find(BinaryImage==1);
    [iRow, iCol, iSli] = ind2sub(size(BinaryImage),iLin);
end

if cropMode == 0
    % segmented and cropped
    SegmentedImage = zeros(size(I));
    SegmentedImage(iLin) = I(iLin);
    SegmentedImage = SegmentedImage(min(iRow):max(iRow), min(iCol):max(iCol), min(iSli):max(iSli));
else
    SegmentedImage = I(min(iRow):max(iRow), min(iCol):max(iCol), min(iSli):max(iSli));
end

% Replace the background with zero pixels
% SegmentedImage = zeros([abs(max(iRow)-min(iRow))+1, abs(max(iCol)-min(iCol))+1, abs(max(iSli)-min(iSli))+1]);
% if(lDisplay), multiWaitbar( 'Cropping', 0 ); end;
% for i=1:1:size(BinaryImage,1)
%     for j=1:1:size(BinaryImage,2)
%         for k=1:1:size(BinaryImage,3)
%             if BinaryImage (i,j,k) == 1;
%                 SegmentedImage(i-min(iRow)+1,j-min(iCol)+1,k-min(iSli)+1) = I(i,j,k);
%             end
%             if(lDisplay), multiWaitbar( 'Cropping', ((i - 1)*size(BinaryImage,2) + (j-1)*size(BinaryImage,3) + k)/(numel(BinaryImage))); end;
%         end
%     end
% end
% if(lDisplay), multiWaitbar( 'Cropping', 'Close'); end;

SegmentedBackground = zeros(size(I));
SegmentedBackground(~BinaryImage) = I(~BinaryImage);

end

