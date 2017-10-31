function [startIndices, endIndices, outputArray] = findRuns(...
    binaryArray, minLength, maxLength, ignoreZeros)
% FINDRUNS finds the start and end indices of runs of consecutive 1's in
% BINARYARRAY that are longer than MINLENGTH (default 0) and shorter
% than MAXLENGTH (default infinity). IGNORENZEROS is an optional integer
% argument; if set, runs of zeros of length IGNOREZEROS or less are
% ignored.
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

szb = size(binaryArray);
binaryArray = binaryArray(:);

% Immediately return empty arrays if input is empty
if isempty(binaryArray)
    startIndices = [];
    endIndices = [];
    return
end

% Turn runs of zeros of length IGNOREZEROS or less into ones
if nargin > 3 && ignoreZeros > 0
    [shortZeroStart, shortZeroEnd] = findRuns(1-binaryArray,0,ignoreZeros);
    for irun = 1:length(shortZeroStart)
        binaryArray(shortZeroStart(irun):shortZeroEnd(irun)) = 1;
    end
end

% Find start/end of runs of 1's by looking for locations where the
% difference between values in BINARYARRAY is +/- 1, also adjusting for the
% start and end of the array.
if binaryArray(1) == 1
    startIndices = [1; find(diff(binaryArray) == 1)+1];
else
    startIndices = find(diff(binaryArray) == 1)+1;
end
if binaryArray(end) == 1
    endIndices = [find(diff(binaryArray) == -1); length(binaryArray)];
else
    endIndices = find(diff(binaryArray) == -1);
end

% Remove runs that are shorter than MINLENGTH
if nargin >= 2
    tooShort = (endIndices - startIndices + 1) < minLength;
    startIndices = startIndices(~tooShort);
    endIndices = endIndices(~tooShort);
end

% Remove runs that are longer than MAXLENGTH
if nargin >= 3 && ~isempty(maxLength) && maxLength > 0
    tooLong = (endIndices - startIndices + 1) > maxLength;
    startIndices = startIndices(~tooLong);
    endIndices = endIndices(~tooLong);
end

% If specified as an output, return the binary array only including valid
% runs
if nargout >= 3
    outputArray = false(size(binaryArray));
    for irun = 1:length(startIndices)
        outputArray(startIndices(irun):endIndices(irun)) = true;
    end
    outputArray = reshape(outputArray, szb);
end