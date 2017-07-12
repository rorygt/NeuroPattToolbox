function [phi, v0, v_direction, fullStruct] = orderParameter(vx, vy)
% Function to calculate order parameter, from equation 1 in Collective
% Motion by Vicsek & Zafeiris.
% VX and VY must be NxMxtime vectors giving the x and y components of the
% vector field. If VX and VY have more than 3 dimensions, extra dimensions
% will be reshaped into the third dimension to calculate statistics.

% Accept a single complex input with real part VX and imaginary VY
if ~isreal(vx) && (~exist('vy', 'var') || isempty(vy))
    vy = imag(vx);
    vx = real(vx);
end


szv = size(vx);
szt = [szv(3:end), 1, 1];

% Check that VX and VY are the same size
if any(size(vx) ~= size(vy))
    error('VX and VY must be the same size!')
elseif ndims(vx) > 3
    % Combine extra dimensions into the third dimension
    vx = vx(:,:,:);
    vy = vy(:,:,:);
end

n = numel(vx(:,:,1));

% Calculate average velocity magnitude
v0full = sqrt(vx.^2 + vy.^2);
v0 = squeeze(mean(nanmean(v0full)));

% Calculate average normalized velocity magnitude
vx_sum = squeeze(sum(nansum(vx)));
vy_sum = squeeze(sum(nansum(vy)));
phi = 1./(n*v0) .* sqrt(vx_sum.^2 + vy_sum.^2);

% Reshape to include more than 3 dimensions if inputs had more
phi = reshape(phi, szt);
v0 = reshape(v0, szt);

% Calculate average velocity direction
if nargout > 2
    meanVect = vx_sum + 1i*vy_sum;
    v_direction = reshape(meanVect./abs(meanVect), szt);
end

% Calculate divergence and curl at each time step
if nargout > 3
    fullStruct.order = phi;
    fullStruct.meanMag = v0;
    fullStruct.meanDir = v_direction;
    
    nt = size(vx, 3);
    vdiv = zeros(nt, 1);
    vcurl = zeros(nt, 1);
    for it = 1:nt
        % Don't include edge rows, as they can skew the statistics
        idiv = divergence(vx(:,:,it), vy(:,:,it));
        idiv = idiv(2:end-1,2:end-1);
        vdiv(it) = nanmedian(idiv(:));
        icurl = curl(vx(:,:,it), vy(:,:,it));
        icurl = icurl(2:end-1,2:end-1);
        vcurl(it) = nanmedian(icurl(:));
    end
    fullStruct.div = reshape(vdiv, szt);
    fullStruct.curl = reshape(vcurl, szt);
end