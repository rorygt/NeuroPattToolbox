function plotPatternLocs(allLocs, pattTypes, realTime, ntrials, nspacebins, smoothScale, useContours)
% Takes the patterns cell array ALLPATTERNLOCS from findAllPatterns.m and
% produces a number of plots

plotTime = true;
if isscalar(realTime)
    realTime = 1:realTime;
    plotTime = false;
end
if ~exist('ntrials', 'var')
    ntrials = 1;
end
if ~exist('smoothScale', 'var')
    smoothScale = 1;
end
if ~exist('useContours', 'var')
    useContours = true;
end

% Combine patterns across all trials
npat = length(pattTypes);
combLocs = cell(1, npat);
noPatts = false(1, npat);
if ~iscell(allLocs{1})
    allLocs = {allLocs};
end

for ipatt = 1:npat
    % Reshape cell array to combine across all trials
    thisPatt = cellfun(@(x) x{ipatt}, allLocs, 'UniformOutput', 0);
    thisPatt = cat(1, thisPatt{:});
    if isempty(thisPatt)
        noPatts(ipatt) = true;
    else
        combLocs{ipatt} = thisPatt;
    end
end

combLocs = combLocs(~noPatts);
pattTypes = pattTypes(~noPatts);
npat = length(pattTypes);

% Plot distributions over time and space
for ipatt = 1:npat
    thisPatt = combLocs{ipatt};
    
    % Take histogram of pattern times
    nActive = histcounts(thisPatt(:,3), 1:length(realTime));
    nActive = nActive / ntrials;
    percActive = 100 * sum(nActive) / length(realTime);
    if plotTime
        subplot(2, npat, ipatt)
        plot(realTime(1:end-1), smooth(nActive, smoothScale))
        axis tight
        xlabel('Time (s)')
        ylabel('# patterns active')
        title(sprintf('%s, %0.1f%% active', pattTypes{ipatt}, percActive))
        subplot(2, npat, npat+ipatt)
    else
        subplot(2, ceil(npat/2), ipatt)
    end
    
    % Take 2D histogram of pattern locations
    thisComplexLoc = thisPatt(:,1) + 1i*thisPatt(:,2);
    if all(abs(abs(thisComplexLoc(:))-1) < 10*eps)
        % Check if pattern is plane wave with direction instead of location
        polarhistogram(angle(thisComplexLoc), 24)
        title('Plane wave propagation directions')
    else
        if useContours
            [counts,xedge,yedge] = histcounts2(thisPatt(:,1), ...
                thisPatt(:,2), nspacebins);
            xs = xedge(1:end-1) + 0.5 * xedge(2) - xedge(1);
            ys = yedge(1:end-1) + 0.5 * yedge(2) - yedge(1);
            contourf(imresize(counts, 2), 10)
            colorbar
            set(gca, 'XTick', [], 'YTick', [])
        else
            histogram2(thisPatt(:,1), thisPatt(:,2), nspacebins, ...
                'FaceColor', 'flat', 'ShowEmptyBins', 'on', ...
                'DisplayStyle', 'tile');
            axis equal
            set(gca, 'XTick', [], 'YTick', [])
            xlabel('X location')
            ylabel('Y location')
            title('Spatial centre locations')
            colorbar
        end
    end
    
    if ~plotTime
        title(sprintf('%s, %0.1f%% active', pattTypes{ipatt}, percActive))
    end
    

    
end

