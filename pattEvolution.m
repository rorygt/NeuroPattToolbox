function [nobs, nexp] = pattEvolution(patterns, nt, nafter, nbefore)
% PATTEVOLUTION finds the number of transitions between pattern types.
% Calculates the number of times in a recording that transitions
% from pattern to pattern were observed (NOBS) and the expected number of
% transitions (NEXP) if patterns occurred randomly in time.
%
% INPUTS: pattern matrix PATTERNS (as output by FINDALLPATTERNS), the total
% number of time steps in the recording NT, and the number of time steps
% after (NAFTER) and before (NBEFORE) the end of one pattern occurence to
% search for transitions.
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

if ~iscell(patterns)
    patterns = {patterns};
end

if ~exist('nbefore', 'var')
    nbefore = 0;
end

ntrials = numel(patterns);
npatts = max(cellfun(@(x) max([0; x(:,1)]), patterns));
nobs = zeros(npatts, npatts, ntrials);
nexp = nobs;

% Iterate over all trial repetitions
for itrial = 1:ntrials
    thisPatts = patterns{itrial}(:,1:3);
    concCounts = zeros(npatts);
    
    % Iterate over every pattern to count consecutive occurences
    for ipatt = 1:size(thisPatts, 1)
        % Find consecutive patterns, given maximum number of time steps
        % before and after to count as consecutive
        thisType = thisPatts(ipatt, 1);
        endStartDiff = thisPatts(:,2) - thisPatts(ipatt,3);
        isConc = endStartDiff <= nafter & endStartDiff >= -nbefore;
        concTypes = thisPatts(isConc, 1);
        concCounts(thisType, :) = concCounts(thisType, :) + ...
            histc(concTypes(:)', 1:npatts);
    end
    nobs(:,:,itrial) = concCounts;
    
    % Calculate expected occurences based on pattern prevalence
    % First find odds of every type of pattern occurring in any overlap
    % window
    totOcc = histc(thisPatts(:,1), 1:npatts);
    windSize = nbefore + nafter + 1;
    diffnexp = totOcc * totOcc' * windSize / nt;
%     % A consecutive occurence of one type of pattern can only happen after
%     % the previous occurence has finished
%     samenexp = totOcc * totOcc' * (nafter-1) / nt;
%     thisExp = diffnexp .* (1-eye(npatts)) + samenexp .* eye(npatts);
    nexp(:,:,itrial) = diffnexp;
    
end
