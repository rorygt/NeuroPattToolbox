% Main script to filter data, calculate optical flow, find patterns, and
% view results
%
% This will probably be changed to a function later, but is currently a
% script for ease of testing. 
%
% Required variables:
%   - LFPs: X x Y x TIME x REPETITION matrix representing some
%     oscillatory recording data. If data is spontaneous or stimulus is
%     not repeated, LFPs should be X x Y x TIME.
%   - Fs: Sampling frequency (in Hz)
%   - dataName: (Optional) String identifying data being used (e.g.
%     experiment and animal ID), will be included in figure titles and
%     saved data

clearvars
close all

% This is a customized function I use to load relevant data, other users
% will need to replace this with some other function to load data
loadDataRory

%% Set parameters
% Most parameters are defined in the SETPARAMS function, this section sets
% secondary parameters only used in this script

% Main set of parameters from dedicated function
params = setParams(Fs);

% Optionally give time limits by which to crop evoked data (leave empty
% for no cropping)
if size(LFPs,4) > 1
    timeCropLims = [-0.5, 1.5];
else
    timeCropLims = [];
end

% Video output parameters
% Flag to save video of signal, velocity fields and patterns
saveVideo = false;
% Video file name (including path). If empty, video will be saved in
% default MATLAB path figure title below as a file name
vidName = [];
% Video frame rate (frames/s)
vidFps = 20;
% Scale by which to resize phase/amplitude maps in video
vidResizeScale = 2;
% Scale by which to increase velocity field vectors in video
vidVectorScale = 1.5;

% Set title for figures and save files
if params.useAmplitude
    typeStr = 'amp';
else
    typeStr = 'phase';
end
if params.useMorlet
    figTitle = sprintf('%s_%iHz_mp%i_a%0.1fb%0.1f_%s', dataName, ...
        params.morletCfreq, params.morletParam, params.opAlpha, ...
        params.opBeta, typeStr);
else
    figTitle = sprintf('%s_%i-%iHz_a%0.1fb%0.1f_%s', dataName, ...
        params.hilbFreqLims(1), params.hilbFreqLims(2), params.opAlpha, ...
        params.opBeta, typeStr);
end

% Some plotting parameters are set within the relevant section below so
% that sections can be run independently for easy plotting. They are
% initially defined here to make them easier to find, but their values here
% are not used.
% Parameters for finding pattern evolutions
nafter = round(0.05*Fs);
nbefore = round(0.01*Fs);
% Parameters for plotting pattern locations
smoothScale = 20;
nspacebins = 2*(9 - 2*params.minEdgeDistance);
% Parameters for plotting velocity field modes
useComplexSVD = false;
vectorScale = 1.5;

%% Process LFPs
% Find any channels with NaN or zero values
nanChans = any(isnan(LFPs(:,:,:)),3);
zeroChans = all(LFPs(:,:,:)==0, 3);
badChannels = find(nanChans | zeroChans);
if ~useMouseCortex
    % Also exclude corner electrodes from Utah MEA LFPs
    badChannels = union([1 10 91 100], badChannels);
end

% Filter LFPs
disp('Filtering waveforms...'); tic
timeDim = 3;
wvcfs = squeeze(morletWaveletTransform(LFPs, Fs, params.morletCfreq, ...
    params.morletParam, timeDim));
toc

% Downsample data
wvcfs = wvcfs(:,:,1:params.downsampleScale:end,:);
Fs = Fs/params.downsampleScale;
if useEvoked
    realTime = (0:(size(wvcfs,timeDim)-1))/Fs-1;
else
    realTime = (1:size(wvcfs,timeDim))/Fs;
end

% Optionally crop data
if ~isempty(timeCropLims)
    cropInds = find(realTime>=timeCropLims(1),1) : ...
        find(realTime>=timeCropLims(2),1);
    wvcfs = wvcfs(:,:,cropInds,:);
    realTime = realTime(cropInds);
end

%% Test the maximum difference between successive time points
if params.useAmplitude
    maxRange = max(abs(wvcfs),[],timeDim) - min(abs(wvcfs),[],timeDim);
    allDiff = abs(diff(abs(wvcfs),1,timeDim)) ./ ...
        repmat(maxRange,1,1,size(wvcfs,3)-1,1);
else
    allDiff = anglesubtract(angle(wvcfs(:,:,2:end,:)), ...
        angle(wvcfs(:,:,1:end-1,:)), 1);
    allDiff = allDiff / pi;
end
maxDiff = prctile(abs(allDiff(:)), 99);
fprintf('99th percentile of the fractional change between time steps is %0.2f.\n', maxDiff)
if maxDiff > 0.1
    disp('Change is >10%, results may be affected by low sampling frequency.')
    if params.downsampleScale > 1
        disp('Consider reducing downsampleScale in setParams.m to increase Fs.')
    end
else
    disp('Change is <10%, sampling frequency is sufficient for data.')
    disp('Consider increasing downsampleScale in setParams.m to improve speed.')
end

clearvars allDiff

%% Compute optical flow
disp('Computing optical flow fields...'); tic

vfs = zeros(size(wvcfs));
vfs = vfs(:,:,1:end-1,:);
meanCSteps = zeros(size(wvcfs,4), 1);

% Calculate velocity fields for every trial, and same average number of
% steps to converge
for itrial = 1:size(wvcfs,4)
    [vx, vy, csteps] = opticalFlow(wvcfs(:,:,:,itrial), badChannels, ...
        params.opAlpha, params.opBeta, ~params.useAmplitude);
    vfs(:,:,:,itrial) = vx + 1i*vy;
    if useMouseCortex
        vfs(:,:,:,itrial) = vfs(:,:,:,itrial) .* ...
            repmat(cortexMask, 1, 1, size(vfs,3));
    end
    
    meanCSteps(itrial) = mean(csteps);
    fprintf('Processed trial %i\n', itrial)
end

toc
fprintf('Optical flow took %0.1f steps on average to converge.\n', mean(meanCSteps))

%% Loop over every trial to find patterns present
disp('Finding all patterns...'); tic;

allPatts = cell(1, size(wvcfs,4));
allLocs = allPatts;

for itrial = 1:size(wvcfs,4)
    thisvf = vfs(:,:,:,itrial);
    [patterns, pattTypes, colNames, pattLocs] = ...
        findAllPatterns(real(thisvf), imag(thisvf), params, ...
        angle(wvcfs(:,:,:,itrial)));
    allPatts{itrial} = patterns;
    allLocs{itrial} = pattLocs;
end

activeArray = makeActivePatternsArray(allPatts, length(pattTypes), size(wvcfs,3));
figure
imagesc((1:size(wvcfs,3))/Fs, 1:length(pattTypes), activeArray)
set(gca, 'YTick', 1:length(pattTypes), 'YTickLabel', pattTypes)

toc

%% Examine evolution between patterns
% Number of time steps before and after a pattern ends to search for other
% patterns
nafter = round(0.05*Fs);
nbefore = round(0.01*Fs);

pattTypeStr = '';
for itype = 1:length(pattTypes)
    pattTypeStr = sprintf('%s%i.%s ', pattTypeStr,itype,pattTypes{itype});
end

[nobs, nexp] = pattEvolution(allPatts, length(realTime), nafter, nbefore);
rateDiff = (nobs - nexp) / length(realTime) * Fs;
%rateDiff = (nobs - nexp);% ./ (nexp);
disp('Observed minus expected pattern transitions/sec')
disp(pattTypeStr)
disp(nanmean(rateDiff,3))
%disp(median(nobs,3) - median(nexp,3))

disp('Wilcoxon sign rank test p-values')
pvals = zeros(size(nobs,1));
for initPatt = 1:size(nobs,1)
    for nextPatt = 1:size(nobs,2)
        thisObs = nobs(initPatt, nextPatt, :);
        thisExp = nexp(initPatt, nextPatt, :);
        [h, p] = ttest(thisObs(:),  thisExp(:));
        pvals(initPatt, nextPatt) = p;
    end
end
disp(pvals)

% symbolArray = zeros(1, size(activeArray,2));
% for ipatt = 1:size(activeArray,1)
%     symbolArray(activeArray(ipatt,:)>0) = ipatt;
% end
% motifs = temporalMotif(symbolArray, [2 3 4 5], 0);

%% Plot pattern locations
figure('Name', figTitle)
% Optionally smooth pattern counts in time
smoothScale = 20;
% Number of bins in x- and y-directions to plot pattern location counts
nspacebins = 2*(9 - 2*params.minEdgeDistance);
plotPatternLocs(allLocs, pattTypes, realTime, size(wvcfs,4), nspacebins, smoothScale, 1);

%% Plot patterns duration and displacement distributions
figure('Name', figTitle)
totPatts = cat(1, allPatts{:});
unqPatts = unique(totPatts(:,1));
maxDur = max(totPatts(:, strcmp(colNames, 'duration')))/Fs;
maxDisp = max(totPatts(:, strcmp(colNames, 'meanDisplacement')))*Fs;
for ipatt = 1:length(unqPatts)
    thisPatt = totPatts(totPatts(:,1)==unqPatts(ipatt), :);
    % Histogram of durations
    subplot(2, length(unqPatts), ipatt)
    thisDurs = thisPatt(:, strcmp(colNames, 'duration'));
    histogram(thisDurs/Fs, linspace(0, maxDur, 10))
    xlabel('Duration (s)')
    ylabel('Counts')
    title(sprintf('%s, mean %0.3g', pattTypes{unqPatts(ipatt)}, mean(thisDurs/Fs)))
    
    % Histogram of displacements
    subplot(2, length(unqPatts), length(unqPatts)+ipatt)
    thisDisp = thisPatt(:, strcmp(colNames, 'meanDisplacement'));
    histogram(thisDisp*Fs, linspace(0, maxDisp, 10))
    xlabel('Displacement (grid/s)')
    ylabel('Counts')
    title(sprintf('Mean %0.3g', mean(thisDisp*Fs)))
end

%% Perform singular value decomposition of velocity fields
% Flag to use complex SVD, which means that spatial modes are free to
% rotate (otherwise they can only be scaled)
useComplexSVD = false;
% Scale by which to increase vector length from MATLAB default
vectorScale = 1.5;

% Open new figure and plot SVD modes
figure('Name', figTitle)
if useEvoked
    plotTime = realTime(1:end-1);
else
    plotTime = [];
end
plotcsvd(vfs, 6, plotTime, useComplexSVD, vectorScale);

%% Optionally save a video file of all data
if saveVideo
    % Set video name
    if isempty(vidName)
        vidName = strrep(figTitle, '.', '');
    end
    
    saveVelocityFieldVideo(wvcfs, vfs, vidName, vidFps, ...
        Fs, vidResizeScale, vidVectorScale, params.useAmplitude)
    
end

 

