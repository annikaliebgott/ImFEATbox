function TransformedData = fEstimateReduction_newData(Samples,X,mappedX)
% This function estimated the reduced feature matrix of new data from the
% original feature matrix and the corresponding calculated reduced feature
% matrix.
%
% ************************************************************************
% Modified for MRI feature extraction by the Department of Diagnostic
% and Interventional Radiology, University Hospital of Tuebingen, Germany
% and the Institute of Signal Processing and System Theory University of
% Stuttgart, Germany. Last modified: September 2018
%
% This implementation is part of ImFEATbox, a toolbox for image feature
% extraction and analysis. Available online at:
% https://github.com/annikaliebgott/ImFEATbox
%
% Contact: annika.liebgott@iss.uni-stuttgart.de
% ************************************************************************
%
% Implementation by: (C) Laurens van der Maaten,
%                        Delft University of Technology


% Remove duplicates from the dataset
X = double(unique(X, 'rows'));

% For all datapoints
TransformedData = zeros(size(Samples, 1), size(mappedX, 2));
bb = sum(X' .* X');
for i=1:size(Samples, 1)
    
    % Get current point
    point = Samples(i,:);
    
    % Find nearest neighbor for current point
    aa = sum(point .* point);
    ab = point * X';
    d = sqrt(repmat(aa', [1 size(bb, 2)]) + repmat(bb, [size(aa, 2) 1]) - 2 * ab);
    [~, ind] = min(d);
    
    % Compute transformation matrix
    L = pinv(X(ind,:) - mean(X(ind,:))) * (mappedX(ind,:) - mean(mappedX(ind,:)));
    
    % Compute coordinates of transformed point
    TransformedData(i,:) = (mean(mappedX(ind,:)) + ((point - mean(X(ind,:))) * L))';
end
end