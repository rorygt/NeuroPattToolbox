function params = setPatternParams(type, Fs)
% Function to set parameters for pattern detection, given sampling
% frequency Fs (in Hz)

%% Filtering parameters
% Optionally downsample after filtering for faster optical flow and pattern
% detection calculations
downsampleScale = 5;
% Morlet transform parameters
% Centre frequency (in Hz)
morletCfreq = 6;
% Morlet wavelet scale parameter: higher values have better frequency
% resolution but lower temporal resolution
morletParam = 5;
% Hilbert transform parameters
% Frequency limits of band-pass filter (in Hz)
hilbFreqLims = [1 4];

%% Optical flow parameters
% Smoothing parameter: higher values will give a smoother velocity field
% (typically 0<OPALPHA<5).
opAlpha = 0.5;
% Non-linearity penalty parameter: close to zero will be highly non-linear,
% large values will be approximately linear, resulting in faster
% computations but possibly less accurate flow fields
opBeta = 1;
% Flag to calculate amplitude velocity fields rather than phase
useAmplitude = false;

%% Pattern parameters
% Minimum order parameter for a plane wave to be detected
planeWaveThreshold = 0.85;
% Maximum velocity field magnitude for synchrony to be detected
synchronyThreshold = 0.85;
% Minimum duration of a pattern for it to be stored (in seconds)
minDurationSecs = 0.02;
% Maximum duration between critical points (or synchrony/plane
% waves) for them to be counted as the same pattern (in seconds)
maxTimeGapSecs = 0.005;
% Maxiumum displacement between critical points between time steps for them
% to be counted as the same pattern (measured in grid spaces)
maxDisplacement = 1;
% Minimum spatial radius for a critical point to occupy for it to be
% counted, quantified by the winding number (measured in grid spaces)
minCritRadius = 2;
% Minimum distance from the edge of the system (in grid spaces)
minEdgeDistance = 2;
% Boolean paramter to combine node and focus type critical points
combineNodeFocus = false;
% Boolean parameter to combine stable and unstable critical points
combineStableUnstable = false;

%% Set parameters in structure based on input type
% Find sampling frequency
if nargin==1 && isnumeric(type)
    Fs = type;
end
% Set all variables by default
if nargin==0 || isempty(type) || isnumeric(type)
    type = 'all';
end

if strcmp(type, 'all') || strcmp(type, 'filt')
    % Filtering parameters
    params.downsampleScale = downsampleScale;
    params.morletCfreq = morletCfreq;
    params.morletParam = morletParam;
    params.hilbFreqLims = hilbFreqLims;
end
if strcmp(type, 'all') || strcmp(type, 'optic')
    % Optical flow parameters
    params.opAlpha = opAlpha;
    params.opBeta = opBeta;
    params.useAmplitude = useAmplitude;
end
if strcmp(type, 'all') || strcmp(type, 'pattern')
    % Pattern detection parameters
    params.planeWaveThreshold = planeWaveThreshold;
    params.synchronyThreshold = synchronyThreshold;
    params.minDurationSecs = minDurationSecs;
    params.maxTimeGapSecs = maxTimeGapSecs;
    params.maxDisplacement = maxDisplacement;
    params.minCritRadius = minCritRadius;
    params.minEdgeDistance = minEdgeDistance;
    params.combineNodeFocus = combineNodeFocus;
    params.combineStableUnstable = combineStableUnstable;
    
    % Set parameters using sampling frequency
    if ~exist('Fs', 'var') || ~isnumeric(Fs) || Fs<=0
        error('Pattern detection parameters can only be set if sampling frequency is input!')
    end
    params.minDuration = max(1, round(params.minDurationSecs * Fs));
    params.maxTimeGap = floor(params.maxTimeGapSecs * Fs);
end