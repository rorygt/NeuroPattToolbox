function [velocityX, velocityY, allConvSteps] = ...
    opticalFlow(video, badChannels, alpha, beta, angleFlag, nStepsDisplay)
% OPTICALFLOW estimates the velocity of pixels between frames in the 
% (x x y x t) matrix VIDEO. 
%
% INPUTS:
%   VIDEO: (X x Y x T) matrix of successive frames of intensity or phase
%       data
%   BADCHANNELS: vector containing the indices of all channels to be
%       ignored in VIDEO.
%   ALPHA: smoothness weighting parameter, larger values of ALPHA will
%       create a smoother optical flow field (typically 0<ALPHA<5).
%   BETA: nonlinear penalty parameter. If BETA is close to zero, optical
%       flow will be computed with highly nonlinear penalties, usually
%       resulting in a more accurate field which is slower to compute.
%   ANGLEFLAG: flag indicating that VIDEO contains angular (phase) data.
%   NSTEPSDISPLAY: optional, integer specifying the number of steps
%       before a message is displayed. If NSTEPSDISPLAY is not defined or
%       is set to 0, no messages will be displayed, otherwise the progress
%       of the calculation will be updated every NSTEPSDISPLAY steps.
%
% Rory Townsend, Aug 2017
% rory.townsend@sydney.edu.au

% Set default inputs
if ~exist('badChannels', 'var')
    badChannels = find(any(isnan(video), 3));
    % Take the gradient of one snapshot, as other channels can have NaN
    % gradient if they are sandwiched between other NaN's
    [gradx, grady] = phasegradient(video(:,:,1), badChannels, 0);
    badChannels = find(isnan(gradx) | isnan(grady));
end
if ~exist('alpha', 'var')
    alpha = 0.2;
end
if ~exist('beta', 'var')
    beta = 1;
end
if ~exist('angleFlag', 'var')
    if ~isreal(video)
        angleFlag = 1;
        disp('OpticalFlow: Assuming that input video contains circular data.')
    else
        angleFlag = 0;
        disp('OpticalFlow: Assuming that input video contains linear data.')
    end
end

% Initialize result structures
nrows = size(video, 1);
ncols = size(video, 2);
nframes = size(video, 3);
ivx = zeros(nrows, ncols);
ivy = ivx;
velocityX = zeros(nrows, ncols, nframes-1);
velocityY = velocityX;
allConvSteps = nan(1, nframes-1);

% If data is not angular, normalize values by scaling by the overall mean
if ~angleFlag
    video = video / nanmean(abs(video(:)));
    % Also convert complex numbers to amplitudes
    if ~isreal(video)
        video = abs(video);
    end
end

% Interpolate over bad electrodes
video = interpolateDeadElectrodes(video, badChannels);

% Initialize temporary variables
prevFrame = [];
frame = video(:,:,1);
nextFrame = video(:,:,2);
if nframes >= 3
    next2Frame = video(:,:,3);
else
    next2Frame = [];
end

% To construct the coefficient matrix A for the system of linear equations
% in OPTICALFLOWSTEP, it is necessary to be able to find the location of
% the surrounding pixels around any given pixel in index form in order to
% estimate spatial derivatives and the laplacian. Set these matrices up now
% to save time inside the loop.
dxMatrix = sparse(nrows*ncols, nrows*ncols);
dyMatrix = dxMatrix;
lapMatrix = dxMatrix;
for irow = 1:nrows
    for icol = 1:ncols
        % Find all surrounding locations in index form
        xSurrLocs = [irow-2, irow-1, irow, irow+1, irow+2] + ...
            nrows * (icol-1) * ones(1,5);
        ySurrLocs = irow * ones(1,5) + ...
            nrows * ([icol-2, icol-1, icol, icol+1, icol+2] - 1);
        
        % Initialize weights of surrounding pixels for laplacian
        lapWeight = [1 1 1 1];
        
        % Deal with cases where the row is close to the edge
        if irow == 1
            % Forward difference
            dxLocs = xSurrLocs([3,4]);
            dxWeight = [-1 1];
            lapWeight = [0 2 1 1];
        elseif irow == nrows
            % Backward difference
            dxLocs = xSurrLocs([2,3]);
            dxWeight = [-1 1];
            lapWeight = [2 0 1 1];
        elseif irow == 2 || irow == nrows-1
            % Centred difference
            dxLocs = xSurrLocs([2,4]);
            dxWeight = [-0.5 0.5];
        else
            % 5-point stencil
            dxLocs = xSurrLocs([1,2,4,5]);
            dxWeight = 1/12 * [1 -8 8 -1];
        end
            
        % Deal with cases where the column is close to the edge
        if icol == 1
            % Forward difference
            dyLocs = ySurrLocs([3,4]);
            dyWeight = [-1 1];
            lapWeight = lapWeight + [0 0 -1 1];
        elseif icol == ncols
            % Backward difference
            dyLocs = ySurrLocs([2,3]);
            dyWeight = [-1 1];
            lapWeight = lapWeight + [0 0 1 -1];
        elseif icol == 2 || icol == ncols-1
            % Centred difference
            dyLocs = ySurrLocs([2,4]);
            dyWeight = [-0.5 0.5];
        else
            % 5-point stencil
            dyLocs = ySurrLocs([1,2,4,5]);
            dyWeight = 1/12 * [1 -8 8 -1];
        end  
            
        % Choose locations and weights for the laplacian
        lapLocs = [xSurrLocs([2,4]), ySurrLocs([2,4])];
        lapLocs = lapLocs(lapWeight~=0);
        lapWeight = lapWeight(lapWeight~=0);
        
        % Save surrounding locations in matrices
        thisRow = sub2ind([nrows, ncols], irow, icol);
        dxMatrix(thisRow, dxLocs) = dxWeight;
        dyMatrix(thisRow, dyLocs) = dyWeight;
        lapMatrix(thisRow, lapLocs) = lapWeight;
    end
end
lapMatrix = lapMatrix + diag(-4*ones(nrows*ncols, 1));
surroundLocs.dx = dxMatrix;
surroundLocs.dy = dyMatrix;
surroundLocs.laplacian = lapMatrix;

% Loop over all time steps
for it = 1 : ( size(video, 3) - 1 )
    % Calculate optical flow
    [ivx, ivy, convSteps] = opticalFlowStep(frame, nextFrame, ...
        badChannels, surroundLocs, alpha, ...
        beta, 0, ivx, ivy, prevFrame, next2Frame, angleFlag);
    
    if convSteps == 1000
        disp('HERE')
    end
    
    % Store results
    allConvSteps(it) = convSteps;
    velocityX(:,:,it) = ivx;
    velocityY(:,:,it) = ivy;
    
    % Display the current step every MSTEPSDISPLAY steps
    if exist('nStepsDisplay', 'var') && nStepsDisplay > 0
        if mod(it, nStepsDisplay) == 0
            fprintf('Calculating velocity, step %d\n', it)
        end
    end
    
    % Next set of frames
    prevFrame = frame;
    frame = nextFrame;
    nextFrame = next2Frame;
    if it+3 <= size(video, 3)
        next2Frame = video(:, :, it+3);
    else
        next2Frame = [];
    end
    
end

% Display warning for number of steps that didn't converge
if sum(allConvSteps==1000) > 0
    fprintf('Warning! %i of %i time steps did not converge.\n', ...
        sum(allConvSteps==1000), length(allConvSteps))
end
