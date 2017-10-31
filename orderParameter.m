function [order, meanMag, meanDir, fullStruct] = orderParameter(vx, vy)
% ORDERPARAMETER calculates order parameters of a vector field.
%
% INPUTS: VX and VY are NxMxtime vectors giving the x and y components of
%   the vector field. If VX and VY have more than 3 dimensions, extra
%   dimensions will be reshaped into the third dimension to calculate
%   statistics. Alternately, input can be a single vector where the real
%   component is VX and imaginary component is VY.
%
% OUTPUTS: 
% ORDER is the average normalized velocity magnitude (used as an order
%   parameter in T. Vicsek & A. Zafeiris (2012), Collective Motion, Physics
%   Reports 517:71-140.
% MEANMAG is the average velocity magnitude.
% MEANDIR is the average normalized velocity direction.
% FULLSTRUCT is a structure containing all previous outputs, as well as the
%   average divergence and curl of the vector field.
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

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
meanMag = squeeze(mean(nanmean(v0full)));

% Calculate average normalized velocity magnitude
vx_sum = squeeze(sum(nansum(vx)));
vy_sum = squeeze(sum(nansum(vy)));
order = 1./(n*meanMag) .* sqrt(vx_sum.^2 + vy_sum.^2);

% Reshape to include more than 3 dimensions if inputs had more
order = reshape(order, szt);
meanMag = reshape(meanMag, szt);

% Calculate average velocity direction
if nargout > 2
    meanVect = vx_sum + 1i*vy_sum;
    meanDir = reshape(meanVect./abs(meanVect), szt);
end

% Calculate divergence and curl at each time step
if nargout > 3
    fullStruct.order = order;
    fullStruct.meanMag = meanMag;
    fullStruct.meanDir = meanDir;
    
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