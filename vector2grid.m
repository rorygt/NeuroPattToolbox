function grid = vector2grid(vector)
% Converts a 1x100, 100x1 or 100xt vector into a 10x10 or 10x10xt grid,
% rotated correctly. channelDim optionally specifies the dimension
% specifying channel number, which must have length 100.
%
% NOTE: THIS TRANSFORMATION CHANGES THE INDEXING OF THE DATA
% A row vector is stored as [1 2 3 4 ... 99 100]. However, in a 10x10 grid
% format, Paul's group stores the configuration of electrodes as:
%       [100  ......  10 ]
%       | .            . |
%       | .            . |
%       | .            2 |
%       [91   ......   1 ]
% This is 180 degrees rotated from the default MATLAB reshaping. The
% relationship between the index of a value is a row or column vector IV
% and the index of the same value in this rotated matrix IM is:
%    IV = 101 - IM
sv = size(vector);
nrow = sv(1);
ncol = sv(2);
ngrid = 10;

% Check which dimension has length 100 (preferentially choosing the first
% dimension)
if nrow == 100
    t = sv(end);
elseif ncol == 100
    t = nrow;
    vector = vector';
elseif nrow == 10 && ncol == 10
    return
elseif sqrt(nrow) == floor(sqrt(nrow))
    % NEW FUNCTIONALITY: If none of these cases apply, use the first
    % dimension as indicating channel number, as long at is a square
    t = ncol;
    ngrid = sqrt(nrow);
else
    error('vectorLength', 'Invalid input vector size!')
end


% Reshape to 10x10 grid and rotate
if t==1
    grid = rot90(reshape(vector,ngrid,ngrid), 2);
else
    % Reshape to 10x10 grid at every time step
    grid = zeros([ngrid,ngrid,sv(2:end)]);
    for it = 1:prod(sv(2:end))
        grid(:,:,it) = rot90(reshape(vector(:,it),ngrid,ngrid), 2);
    end
end