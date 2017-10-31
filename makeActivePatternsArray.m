function activeArray = makeActivePatternsArray(patterns, npat, nt)
%MAKEACTIVEPATTERNSARRAY converts pattern matrix to activity array.
%   This function takes the PATTERNS matrix output by FINDALLPATTERNS
%   giving the start and end of each pattern, and converts to a PxT matrix
%   giving the number of patterns of each time present at each time step,
%   where P is the number of pattern types (given by optional parameter
%   NPAT), and T is the total number of time steps (given by optional
%   parameter NT).
%
%   If PATTERNS is a cell array (representing multiple stimulus trials),
%   ACTIVEARRAY will be summed across all cells.
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

% Convert non-cell array to cell array for consistency
if ~iscell(patterns)
    patterns = {patterns};
end

% Set default number of patterns as maximum in array
if ~exist('npat', 'var') || isempty(npat) || npat==0
    npat = max(cellfun(@(x) max(x(:,1)), patterns));
end

% Set default number of time steps as maximum in array
if ~exist('nt', 'var') || isempty(nt) || nt==0
    nt = max(cellfun(@(x) max(x(:,3)), patterns));
end

% Loop over each cell
activeArray = zeros(npat, nt);
for icell = 1:length(patterns)
    thisCell = patterns{icell};
    % Loop over each pattern
    for ipatt = 1:size(thisCell,1)
        thisPatt = thisCell(ipatt,:);
        activeArray(thisPatt(1), thisPatt(2):thisPatt(3)) = ...
            activeArray(thisPatt(1), thisPatt(2):thisPatt(3)) + 1;
    end
end

end

