% Main function to load raw data, filter, calculate optical flow and find
% patterns

clearvars
close all

%% Data file parameters
useMouseCortex = false;
mouseResize = 0.2;
useEvoked = true;

if useEvoked
    fileInd = 3;
    % Set path to stimulus-evoked files
    recordingsEv = {'MA026-14','MA027-8', 'MA026-44', 'MA027-7', ...
        'MY144-111', 'MY147-31', 'MA026-46'};
    if isunix
        dataLoc = './Processed_data/Evoked/';
    else
        dataLoc = 'D:\Evoked';
    end
    dataName = sprintf('evoked_%s.mat', recordingsEv{fileInd});
    stimDir = 2;
else
    fileInd = 2;
    % Set path to spontaneous files
    animalsSp = {'my144a', 'MY147', 'MA026a'};
    fileNumsSp = [101, 53, 5];
    filePartSp = '';
    if isunix
        dataLoc = './Processed_data/';
    else
        error('Data missing!')
    end
    dataName = sprintf('LFPsHilbertSpikes_%s%s-%i.mat', ...
        animalsSp{fileInd}, filePartSp, fileNumsSp(fileInd));
end

if useMouseCortex
    dataName = 'mouseCortex';
end

%% Set parameters
% Main set of parameters from function
params = setPatternParams('all', Fs);

% Optionally include time limits by which to crop evoked data (leave empty
% for no cropping)
if useEvoked && ~useMouseCortex
    timeCropLims = [-0.5, 1.5];
else
    timeCropLims = [];
end

% Video output parameters
% Flag to save video of signal, velocity fields and patterns
saveVideo = false;
% Video file name
vidName = [];
% Video frame rate
vidFrameRate = 20;

%% Load LFPs
if params.useAmplitude
    typeStr = 'amp';
else
    typeStr = 'phase';
end

figTitle = sprintf('%s_%iHz_mp%i_a%0.1fb%0.1f_%s', dataName, ...
    params.morletCfreq, params.morletParam, params.opAlpha, ...
    params.opBeta, typeStr);

% LFPs should be of the form CHANNEL x TIME x TRIAL (only 1 trial for
% spontaneous recordings)
if useMouseCortex
    % Load mouse data file
    load('./Exp011_Fluo_001_part_2_sequenceData_ratioFiltered.mat')
    Fs = 50;
    data = rot90(data(20:115, 81:144, :), 2);
    firstFrame = imresize(data(:,:,1), mouseResize);
    LFPs = zeros([size(firstFrame), size(data, 3)]);
    cortexMask = imresize(rot90(meta.traceMask(20:115, 81:144,:),2), ...
        mouseResize);
    for itime = 1:size(data,3)
        LFPs(:,:,itime) = imresize(data(:,:,itime), mouseResize);
    end
    clearvars data

else
    fprintf('Loading file %s\n', dataName); tic
    if useEvoked
        % Load stimulus-evoked LFP data file
        load(fullfile(dataLoc, dataName), 'allLFPs', 'Fs')
        LFPs = allLFPs{stimDir};
        clearvars allLFPs
        LFPs = permute(LFPs, [1 3 2]);
    else
        % Load spontaneous LFP data file
        load(fullfile(dataLoc, dataName), 'LFPs', 'Fs')
    end
    
    % Convert to X-Y grid of electrodes
    LFPs = vector2grid(LFPs);
end
toc


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

params = setPatternParams([],Fs);

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
nafter = round(0.05*Fs);
nbefore = round(0.01*Fs);

pattTypeStr = '';
for itype = 1:length(pattTypes)
    pattTypeStr = sprintf('%s%i.%s ', pattTypeStr,itype,pattTypes{itype});
end

[nobs, nexp] = pattEvolution(allPatts, length(realTime), nafter, nbefore);
%rateDiff = (nobs - nexp) / length(realTime) * Fs;
rateDiff = (nobs - nexp);% ./ (nexp);
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
smoothScale = 20;
nspacebins = 2*(9 - 2*params.minEdgeDist);
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

%% Perform SVD of vector fields
useComplexSVD = false;
vectorScale = 1.5;
figure('Name', figTitle)
if useEvoked
    plotTime = realTime(1:end-1);
else
    plotTime = [];
end
plotcsvd(vfs, 6, plotTime, useComplexSVD, vectorScale);

%% Optionally save a video file of all data
if useMouseCortex && saveVideo
    combLocs = allLocs{1};
    combLocs = sortrows(cat(1, combLocs{2:end}), 3);
    
    fig = figure('Name', figTitle);
    
    if isempty(vidName)
        vidName = strrep(figTitle, '.', '');
    end
    
    % Set title and open video file
    vidTitle = strcat(vidName, datestr(now,'mmdd'), '_', ...
        datestr(now,'HHMM'), '.avi');
    vidObj = VideoWriter(vidTitle);
    vidObj.FrameRate = vidFrameRate;
    open(vidObj);
    
    
    for itime = 1:size(vfs,3)
        showGridVectorPatterns(angle(wvcfs(:,:,itime)).*cortexMask ...
            + (1-cortexMask)*1.5 ... % Just to change the mask color
            , vfs(:,:,itime), combLocs(combLocs(:,3)==itime,:), [], 2, 1)
        title(sprintf('%0.2f s', itime/Fs))
        writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
    end
    
    close(vidObj);
end



