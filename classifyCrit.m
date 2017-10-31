function [rowcoords, colcoords, patType, jacobians] = classifyCrit(vx, vy, edgeDistance)
% CLASSIFYCRIT finds and classifies the critical points in the vector field
%   defined by VX and VY. Outputs 2 1xN vectors (where N is the number of
%   critical points detected) expressing the row and  column coordinates of
%   each point, a 1xN vector of the type of critical point present, and
%   optionally a 2x2xN matrix expressing the estimated Jacobian at each
%   point. Critical points fewer than EDGEDISTANCE points from the edge of
%   the array will be discarded (defaults to zero).
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

if ~exist('edgeDistance', 'var')
    edgeDistance = 0;
end

% Find critical points
[rowcoords, colcoords] = criticalpointsbilinear(vx, vy, edgeDistance);

jacobians = zeros(2,2,length(rowcoords));
patType = cell(size(rowcoords));

%     vamp = sqrt(vx.^2 + vy.^2);
%     vx = vx ./ vamp;
%     vy = vy ./ vamp;

if ~isempty(rowcoords)
    [vxx, vxy] = gradient(vx);
    [vyx, vyy] = gradient(vy);
end

for ic = 1:length(rowcoords)
    % Find partial derivatives of the 4 corners that the critical point
    % resides in
    ix = rowcoords(ic);
    iy = colcoords(ic);
    %     corners = sub2ind(size(vx), ...
    %         [floor(ix) floor(ix); ceil(ix) ceil(ix)], ...
    %         [floor(iy) ceil(iy); floor(iy) ceil(iy)]);
    corners = sub2ind(size(vx), ...
        [floor(ix) ceil(ix); floor(ix) ceil(ix)], ...
        [ceil(iy) ceil(iy); floor(iy) floor(iy)]);
    
    % Estimate Jacobian at the critical point through bilinear
    % interpolation
    
    % Use corner points in pre-calculated gradients
    ixdec = ix - floor(ix);
    iydec = iy - floor(iy);
    dxx = [1-ixdec ixdec] * vxx(corners) * [1-iydec iydec]';
    dxy = [1-ixdec ixdec] * vxy(corners) * [1-iydec iydec]';
    dyx = [1-ixdec ixdec] * vyx(corners) * [1-iydec iydec]';
    dyy = [1-ixdec ixdec] * vyy(corners) * [1-iydec iydec]';
    ijac = [dxx dxy; dyx dyy];
    jacobians(:,:,ic) = ijac;
    
    %     % Calculate gradient only for corners of critical points
    %     cornersGradx = singleanglegradientnan(vx, corners, 0);
    %     cornersGrady = singleanglegradientnan(vy, corners, 0);
    %     ixdec = ix - floor(ix);
    %     iydec = iy - floor(iy);
    %     dxx = [1-ixdec ixdec] * cornersGradx(:,:,1) * [1-iydec iydec]';
    %     dxy = [1-ixdec ixdec] * cornersGradx(:,:,2) * [1-iydec iydec]';
    %     dyx = [1-ixdec ixdec] * cornersGrady(:,:,1) * [1-iydec iydec]';
    %     dyy = [1-ixdec ixdec] * cornersGrady(:,:,2) * [1-iydec iydec]';
    %     jacobians(:,:,ic) = [dxx dxy; dyx dyy];
    
    %     % Faster but less accurate: Estimate Jacobian for the cell from only
    %     % the corners
    %     dxx = mean(vx(corners(:,2)) - vx(corners(:,1)));
    %     dxy = mean(vx(corners(2,:)) - vx(corners(1,:)));
    %     dyx = mean(vy(corners(:,2)) - vy(corners(:,1)));
    %     dyy = mean(vy(corners(2,:)) - vy(corners(1,:)));
    %     jacobians(:,:,ic) = [dxx dxy; dyx dyy];
    
    % Classify critical point by its Jacobian
    if det(ijac) < 0
        % Saddle point
        itype = 'saddle';
    elseif trace(ijac)^2 > 4*det(ijac)
        if trace(ijac) < 0
            % Stable node
            itype = 'stableNode';
        else
            % Unstable node
            itype = 'unstableNode';
        end
    else
        if trace(ijac) < 0
            % Stable focus
            itype = 'stableFocus';
        else
            % Unstable focus
            itype = 'unstableFocus';
        end
    end
    patType{ic} = itype;

end
    
end