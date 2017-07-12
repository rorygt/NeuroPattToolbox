function [U, S, V] = plotcsvd(vf, nmodeplot, realTime)
% Performs a complex singular vector decomposition of complex velocity
% field VF, which must be a 3D or 4D vector containing
% [row,column,time,trial]. Plots the top NMODEPLOT spatial modes, energies
% and time courses (default 5). REALTIME is an optional vector containing
% the actual time points in VF in seconds. If realTime==0, only spatial
% modes will be plotted

if ~exist('nmodeplot', 'var')
    nmodeplot = 5;
end

if ~exist('realTime', 'var')
    realTime = 1:size(vf,3);
end

if isempty(realTime) || (~isvector(realTime) && realTime == 0)
    onlySpatial = true;
else
    onlySpatial = false;
end

[nr, nc, nt] = size(vf);

% Perform SVD with complex velocity field
svdmat = reshape(vf, nr*nc, [])';

% Separate x and y components into separate variables for SVD
%svdmat = cat(2, real(svdmat), imag(svdmat));

[U,S,V] = svd(svdmat, 0);

%     % Make SVD with specified number of modes
%     nmodescalc = 20;
%     [U,S,V] = svds(svdmat, nmodescalc);

if size(V,1) > nr*nc
    V = V(1:nr*nc,:) + 1i*V(nr*nc+(1:nr*nc),:);
end

% Calculate total energy then crop to specified number of modes
eigVals = diag(S);
totVar = sum(eigVals.^2);
prctVar = 100 * eigVals.^2 / totVar;
U = -U(:,1:nmodeplot);
S = S(1:nmodeplot, 1:nmodeplot);
V = -V(:,1:nmodeplot);

% Find average trajectory across all trials
reUav = zeros(nt, nmodeplot);
imUav = reUav;
absUav = imUav;
for it = 1:nt
    Uav = U(it + (0:nt:size(U,1)-nt), :);
    reUav(it, :) = mean(real(Uav),1);
    imUav(it, :) = mean(imag(Uav),1);
    absUav(it, :) = mean(abs(Uav),1);
end
UavLims = [min([reUav(:); imUav(:); absUav(:)]), ...
    max([reUav(:); imUav(:); absUav(:)])];

% Plot spatial modes containing most energy
for imode=1:nmodeplot
    % Plot spatial mode
    if onlySpatial
        subplot(2, ceil(nmodeplot/2), imode)
    else
        subplot(2,nmodeplot,imode)
    end
    thisMode = reshape(V(:,imode), nr, nc);
    quiver(real(thisMode), imag(thisMode));
    axis equal
    axis(0.5+[0 nr 0 nc])
    title(sprintf('Var = %0.1f%%', prctVar(imode)))
    %title(sprintf('Var = %0.1f%%, %0.2g%%', prctVar(imode), ...
    %    sum(prctVar(1:imode))))
    
    % Plot time course of spatial mode
    if ~onlySpatial
        subplot(2,nmodeplot,nmodeplot+imode)
        plot(realTime, reUav(:,imode), realTime, imUav(:,imode), ...
            realTime, absUav(:,imode))
        ylim(UavLims)
        xlim([realTime(1), realTime(end)])
        xlabel('Time')
    end
end
suptitle('Top SVD modes')