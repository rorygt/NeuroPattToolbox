% Script to load neural recording for processing in the toolbox. This
% script is customised for use specifically for Rory Townsend, future users
% will need to write appropriate methods for their own data format.

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

%% Load data
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