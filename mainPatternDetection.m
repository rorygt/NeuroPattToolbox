% Main function to load raw data, filter, calculate optical flow and find
% patterns

clearvars
close all

%% Data file parameters
fileInd = 1;
useEvoked = true;

if useEvoked
    % Set path to stimulus-evoked files
    recordingsEv = {'MA026-14','MA027-8', 'MA026-44', 'MA027-7', ...
        'MY144-111', 'MY147-31', 'MA026-46'};
    dataLoc = './Processed_data/Evoked/';
    dataLoc = 'E:\Evoked';
    dataName = sprintf('evoked_%s.mat', recordingsEv{fileInd});
    stimDir = 1;
else
    % Set path to spontaneous files
    animalsSp = {'my144', 'MY147'};
    fileNumsSp = [101, 53];
    filePartSp = '';
    dataLoc = './Processed_data/';
    dataName = sprintf('LFPsHilbertSpikes_%s%s-%i.mat', ...
        animalsSp{fileInd}, filePartSp, fileNumsSp(fileInd));
end

%% Filtering parameters
% Centre frequency
cfreq = 10;
% Morlet wavelet scale parameter: higher values have better frequency
% resolution but lower temporal resolution
morletParam = 6;
% Optionally downsample after filtering for faster optical flow and pattern
% detection calculations
downsampleScale = 5;

%% Optical flow parameters
% Smoothing parameter: higher values will give a smoother velocity field
% (typically 0<OPALPHA<5).
opAlpha = 0.2;
% Non-linearity penalty parameter: close to zero will be highly non-linear,
% large values will be approximately linear, resulting in faster
% computations but possibly less accurate flow fields
opBeta = 1;
% Use flag to calculate amplitude velocity fields rather than phase
useAmplitude = false;

%% Load LFPs
disp('Loading data file...'); tic

% LFPs should be of the form CHANNEL x TIME x TRIAL (only 1 trial for
% spontaneous recordings)
if useEvoked
    load(fullfile(dataLoc, dataName), 'allLFPs', 'Fs')
    LFPs = allLFPs{stimDir};
    clearvars allLFPs
    LFPs = permute(LFPs, [1 3 2]);
else
    load(fullfile(dataLoc, dataName), 'LFPs', 'Fs')
end
toc

% Convert to X-Y grid of electrodes
LFPs = vector2grid(LFPs);

%% Process LFPs

% Find any channels with NaN or zero values
nanChans = any(isnan(LFPs(:,:,:)),3);
zeroChans = sum(LFPs(:,:,:),3) == 0;
badChannels = find(nanChans | zeroChans);
% Also exclude corner electrodes from Utah MEA LFPs
badChannels = union([1 10 91 100], badChannels);

% Filter LFPs
disp('Filtering waveforms...'); tic
wvcfs = squeeze(morletWaveletTransform(LFPs, Fs, cfreq, morletParam, 3));
toc

% Downsample data
wvcfs = wvcfs(:,:,1:downsampleScale:end,:);
Fs = Fs/downsampleScale;
if useEvoked
    realTime = (0:(size(wvcfs,3)-1))/Fs-1;
else
    realTime = size(wvcfs,3);
end

% Test the maximum difference between successive time points
if useAmplitude
    maxRange  = max(abs(wvcfs(:))) - min(abs(wvcfs(:)));
    allDiff = diff(abs(wvcfs),1,2);
else
    maxRange = pi;
    allDiff = anglesubtract(angle(wvcfs(:,2:end,:)), ...
        angle(wvcfs(:,1:end-1,:)), 1);
end
maxDiff = prctile(allDiff(:), 99) / maxRange;
fprintf('99th percentile of the fractional change between time steps is %0.1g.\n', maxDiff)

%% Compute optical flow
disp('Computing optical flow fields...'); tic

vfs = zeros(size(wvcfs));
vfs = vfs(:,:,1:end-1,:);
meanCSteps = zeros(size(wvcfs,4), 1);

% Calculate velocity fields for every trial, and same average number of
% steps to converge
for itrial = 1:size(wvcfs,4)
    [vx, vy, csteps] = opticalFlow(wvcfs(:,:,:,itrial), badChannels, ...
        opAlpha, opBeta, ~useAmplitude);
    vfs(:,:,:,itrial) = vx + 1i*vy;
    meanCSteps(itrial) = mean(csteps);
    fprintf('Processed trial %i\n', itrial)
end

toc
fprintf('Optical flow took %0.1f steps on average to converge.\n', mean(meanCSteps))

%% Loop over every trial to find patterns present
disp('Finding all patterns...'); tic;

params = setPatternParams(Fs);

allPatts = cell(1, size(wvcfs,4));
allLocs = allPatts;

for itrial = 1:size(wvcfs,4)
    thisvf = vfs(:,:,:,itrial);
    [patterns, pattTypes, colNames, pattLocs] = ...
        findAllPatterns(real(thisvf), imag(thisvf), params);
    allPatts{itrial} = patterns;
    allLocs{itrial} = pattLocs;
end

toc

%% Plot pattern locations
figure
smoothScale = 20;
nspacebins = 10;
plotPatternLocs(allLocs, pattTypes, realTime, size(wvcfs,4), nspacebins, smoothScale);


%% Perform SVD of vector fields
figure
plotcsvd(vfs, 3, realTime(1:end-1));







