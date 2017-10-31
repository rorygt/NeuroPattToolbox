function [gfx, gfy, badChannels] = phasegradient(f, badChannels, angleFlag, surroundLocs)
% Calculates the phase gradient of complex matrix F
%
% Rory Townsend, Aug 2017
% rory.townsend@sydney.edu.au

sf = size(f);
f = f(:,:,:);

% Set default inputs
if nargin < 2
    badChannels = find(isnan(sum(f, 3)));
end
if nargin < 3
    angleFlag = true;
end
if nargin < 4
    surroundLocs = [];
end
% For now, don't use SURROUNDLOCS: MATLAB's gradient function is faster and
% more accurate
surroundLocs = [];

% Convert phase data to complex coefficients
if angleFlag && isreal(f)
    f = exp(1i*f);
end

% Convert complex data to amplitudes
if ~angleFlag && ~isreal(f)
    f = abs(f);
end


% Find channels with zero amplitude, as these will cause problems
if angleFlag
    badChannels = union(badChannels, find(any(abs(f)==0,3)));
    % Smooth out bad channels
    f = interpolateDeadElectrodes(f, badChannels);
end

% Iteratively calculate gradient
gfx = zeros(size(f));
gfy = gfx;
for it = 1:size(f, 3)
    if angleFlag
        % Use modified MATLAB inbuilt function
        [igfx, igfy] = phasegradientMATLAB(angle(f(:,:,it)), angleFlag);
        gfx(:,:,it) = -igfx;
        gfy(:,:,it) = -igfy;
    else
        % Regular gradient
        [igfx, igfy] = normalGrad(f(:,:,it), surroundLocs);
        gfx(:,:,it) = igfx;
        gfy(:,:,it) = igfy;
    end
end

% Reshape gradients to be the same size as f
gfx = reshape(gfx, sf);
gfy = reshape(gfy, sf);

end

function [frow, fcol] = normalGrad(f, surroundLocs)
% Choose method to calculate regular, non-circular gradient
if ~isempty(surroundLocs)
    % If optional input SURROUNDLOCS is given, use the more accurate stencil
    % provided to calculate gradient (but may be slower)
    frow = reshape(surroundLocs.dy * f(:), size(f));
    fcol = reshape(surroundLocs.dx * f(:), size(f));
    
else
    % Otherwise just use MATLAB's built-in gradient function
    [frow, fcol] = gradient(f);
end

end