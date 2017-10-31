function motifs = temporalMotif(symbolSeries, lengths, weightFlag, removeCycleFlag)
% TEMPORALMOTIF finds temporal motifs in the array SYMBOLSERIES containing
% integers 0 to n representing different behaviours in a time series.
%
% Method taken from "The Motif Tracking Algorithm. Wilson, Birkin and
% Aickelin, 2007".


% Only look at changes to different behaviours
symbolSeries = symbolSeries([true diff(symbolSeries)~=0]);

% Find unique behaviours
symbols = unique(symbolSeries);
n = length(symbols);

% Count number of occurences of each symbol
symCounts = histc(symbolSeries, symbols);
symFract = symCounts / sum(symCounts);

motifs = cell(n^max(lengths), 1);
mindex = 1;
minOcc = 1;

if ~exist('weightFlag', 'var')
    weightFlag = false;
end

if ~exist('removeCycleFlag', 'var')
    removeCycleFlag = false;
end

for ilength = lengths
    % Generate symbol stage matrix: all subsequences of length ilength
    subsequences = findSubsequences(symbolSeries, ilength);
    [uniqueSubs, locs] = unique(subsequences, 'rows');
    
    % Iteratively remove the first location of unique subsequences to weed
    % out all that occur only a few times
    for remove = 1:(minOcc-1)
        subsequences(locs,:) = [];
        [uniqueSubs, locs] = unique(subsequences, 'rows');
    end
    
    % Store motifs if they occur more than once
    for isub = 1:size(uniqueSubs, 1)
        includedSymbols = unique(uniqueSubs(isub,:));
        if length(includedSymbols) > 1 && ~ismember(0, includedSymbols)
            
            % Count occurences of subsequences that are comprised of more
            % than one symbol
            icount = ...
                sum(ismember(subsequences, uniqueSubs(isub,:), 'rows'));
            
            % Adjust for removed occurences
            icount = icount + minOcc - 1;
            
            % Weight by relative occurance of each pattern
            newMotif = uniqueSubs(isub,:);
            if weightFlag
                fracExpect = prod(symFract(newMotif)) / ...
                    prod(1-symFract(newMotif(2:end)));
                numExpect = fracExpect * length(symbolSeries);
                icount = icount / numExpect;
            end
            

            [motifs, mindex] = ...
                   addToMotifList(newMotif, icount, motifs, mindex);
            %end
        end
        
    end
end
motifs = motifs(1:(mindex-1), :);

% Sort motifs by frequency
[~, order] = sort(cellfun(@(x) x(1), motifs), 'descend');
motifs = motifs(order);

% Adds the frequency to the first column of motifs for easy reading and
% remove motifs that involve only 2 pattern types
deleteMotifs = zeros(1, length(motifs));
for im = 1:length(motifs)
    thisMotif = motifs{im}(2:end);
    if length(unique(thisMotif)) < 3
       deleteMotifs(im) = 1;
       % Optionally remove motifs that cycle back to the same pattern with
       % only one pattern in between
     elseif removeCycleFlag && ( any(diff(thisMotif(1:2:end))==0) ||...
             any(diff(thisMotif(2:2:end))==0) )
         deleteMotifs(im) = 1;
    end
end
motifs = motifs(~deleteMotifs);


%% Subfunctions
function subList = findSubsequences(S, sLength)
    % Finds all subsequences of length SLENGTH in array S. SUBLIST is a
    % LISTLENGTH x SLENGTH matrix where each row contains a subsequence at
    % consecutive time steps.
    
    listLength = length(S) - sLength + 1;
    
    subList = zeros(listLength, sLength);
    
    for ii = 1:listLength
        subList(ii,:) = S(ii:(ii+sLength-1));
    end
end

function [motifs, mindex] = addToMotifList(newMotif, count, motifs, mindex)
    % Adds vector NEWMOTIF and occurence count COUNT to MOTIFS cell array,
    % removing existing entries of MOTIFS that are subsequences of
    % NEWMOTIF. MINDEX gives the index of the first non-empty element of
    % MOTIFS.
    
    deleteCells = zeros(size(motifs));
    
    for ii = 1:mindex
        if ~isempty(findstr(motifs{ii}(2:end), newMotif)) && ...
                motifs{ii}(1) <= count
            deleteCells(ii) = 1;
        end
    end
    
    motifs{mindex} = [count newMotif];
    
    motifs = motifs(~deleteCells);
    
    mindex = mindex - sum(deleteCells) + 1;
      
end

end