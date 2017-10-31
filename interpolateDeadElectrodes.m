function signalGrid = interpolateDeadElectrodes(signalGrid, deadElectrodes, angleFlag)
% Interpolate dead electrodes to be the median value of surrounding nodes.
% Input signalGrid can be a 2D matrix with dimension [X Y] or a 3D matrix
% with dimension [X Y T], where t gives time. badElectrodes is
% an array containing the indices of the dead electrodes for the current
% data set.
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

sg = size(signalGrid);
signalGrid = signalGrid(:,:,:);

[nr, nc, nt] = size(signalGrid);
signalGrid = reshape(signalGrid, [], nt);

% If no dead electrodes are supplied, interpolate over corners and any
% channels with NaNs
if nargin == 1
    nanChans = find(any(isnan(signalGrid), 3));
    corners = sub2ind([nr nc], [1 nr 1 nr], [1 1 nc nc]);
    deadElectrodes = union(corners, nanChans);
end
    
ii = 1;
ndead = length(deadElectrodes(:));

deadElectrodes = deadElectrodes(:)';
while ii <= length(deadElectrodes(:))
    id = deadElectrodes(ii);
    ii = ii+1;
    
    % Skip if this electrode is invalid
    if round(id)~=id || id<1 || id>nr*nc
        continue
    end
    
    % Define surrounding nodes in the form [down up left right]
    surrounding = [id+1 id-1 id+nr id-nr];
    
    % Remove nodes that wrap around the edges
    if mod(id,nr) == 0
        surrounding(1) = [];
    elseif mod(id,nr) == 1
        surrounding(2) = [];
    end

    % Remove out of bounds nodes
    outOfBounds = surrounding<1 | surrounding>nr*nc;
    surrounding(outOfBounds) = [];
    
    % Remove surrounding nodes that are corners or other bad nodes unless
    % these are the only available nodes
    badSurround = ismember(surrounding,deadElectrodes);
    if ii<=ndead || sum(badSurround) < length(surrounding)
        surrounding(badSurround) = [];
    end
    
    % If there are no valid surroudning electrodes, push it to the back of
    % the queue
    if isempty(surrounding)
        deadElectrodes = [deadElectrodes, id];
        continue
    end
    
    % Set to median value of surrounding electrodes
    if nargin < 3 || angleFlag ~= 1
        signalGrid(id, :) = nanmedian(signalGrid(surrounding, :), 1);
    else
        signalGrid(id, :) = angle(nanmean(exp(1i*signalGrid(surrounding, :))));
    end
end

% Reshape to original dimensions
signalGrid = reshape(signalGrid, sg);

end