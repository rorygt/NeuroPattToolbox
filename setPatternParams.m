function params = setPatternParams(Fs)
% Function to set parameters for pattern detection, given sampling
% frequency Fs (in Hz)

% Minimum order parameter for a plane wave to be detected
planeWaveThreshold = 0.85;
% Maximum velocity field magnitude for synchrony to be detected
synchronyThreshold = 0;
% Minimum duration of a pattern for it to be stored (in seconds)
minDurationSecs = 0.02;
% Maximum duration between critical points (or synchrony/plane
% waves) for them to be counted as the same pattern (in seconds)
maxTimeGapSecs = 0.002;
% Maxiumum displacement between critical points between time steps for them
% to be counted as the same pattern (measured in grid spaces)
maxDisplacement = 1;
% Minimum spatial radius for a critical point to occupy for it to be
% counted, quantified by the winding number (measured in grid spaces)
minCritRadius = 1;
% Minimum distance from the edge of the system (in grid spaces)
minEdgeDistance = 2;
% Boolean paramter to combine node and focus type critical points
combineNodeFocus = true;

% Set up parameters structure
params.planeWaveThreshold = planeWaveThreshold;
params.synchronyThreshold = synchronyThreshold;
params.maxDisplacement = maxDisplacement;
params.minCritRadius = minCritRadius;
params.minEdgeDist = minEdgeDistance;
params.combineNodeFocus = combineNodeFocus;

% Set parameters using sampling frequency
params.minDuration = round(minDurationSecs * Fs);
params.maxTimeGap = round(maxTimeGapSecs * Fs);