function saveVelocityFieldVideo(wvcfs, vfs, vidName, vidFps, ...
    Fs, resizeScale, vfScale, useAmplitude)
% Function to write amplitude or phase maps and velocity fields to a
% video file
% INPUTS:
%   - wvcfs, XxYxT matrix of complex wavelet or Hilbert coefficients
%   - vfs, XxYxT matrix of velocity field vectors (as complex values)
%   - vidName, the title for the video file (this function will
%       automatically append the current data and time, and the .avi
%       extension)
%   - vidFps, the frame rate (per second) of the video file (default 20).
%   - Fs, optional, the sampling frequency. Used to display the time step
%       at each video frame (default 1).
%   - resizeScale, optional, the scale by which to spatially interpolate
%       video frames at each time step (must be >=1, default 1).
%   - vfScale, optional, the scale by which to extend or shrink all
%       velocity field vectors (default 1).
%   - useAmplitude, optional, boolean flag indicating video should plot
%       amplitude instead of phase (default false)

disp('Writing video to file...'); tic

% Set default inputs
if ~exist('vidFps', 'var') || isempty(vidFps) || vidFps<1
    vidFps = 20;
end
if ~exist('Fs', 'var') || isempty(resizeScale)
    Fs = 1;
    timeUnit = '';
else
    timeUnit = 's';
end
if ~exist('resizeScale', 'var') || isempty(resizeScale) || resizeScale<=0
    resizeScale = 1;
end
if ~exist('vfScale', 'var') || isempty(vfScale) || resizeScale<=0
    vfScale = 1;
end
if ~exist('useAmpitude', 'var') || useAmplitude ~= 1
    useAmplitude = false;
end

% Create video file
fig = figure;
vidTitle = strcat(vidName, datestr(now,'ddmm'), '_', ...
    datestr(now,'HHMM'), '.avi');
vidObj = VideoWriter(vidTitle);
vidObj.FrameRate = vidFps;
open(vidObj);

% Advance time and save video frames
for itime = 1:size(vfs,3)
    
    % Optionally interpolate spatial grid of data
    signalGrid = imresize(wvcfs(:,:,itime), resizeScale);
    
    % Set signal and colors based on whether you are plotting amplitude
    % or phase
    if useAmplitude
        signalGrid = abs(signalGrid);
        colorMapSpec = parula;
        sigLims = [min(abs(wvcfs(:))), max(abs(wvcfs(:)))];
        vfColor = [1 1 1];
    else
        signalGrid = angle(signalGrid);
        colorMapSpec = pmkmp_new;
        sigLims = [-pi pi];
        vfColor = [0 0 0];
    end
    
    % Plot signal grid
    imagesc(signalGrid, sigLims)
    colormap(gca, colorMapSpec)
    
    % Plot velocity field
    hold on
    vf = vfs(:,:,itime) * vfScale;
    quiver(linspace(1, size(signalGrid,2), size(vf,2)), ...
        linspace(1, size(signalGrid,1), size(vf,1)), ...
        real(vf), imag(vf), 0, 'Color', vfColor)
    set(gca,'YDir','reverse');
    hold off
    
    % Update title with current time
    title(sprintf('%g %s', itime/Fs, timeUnit))
    
    % Write to video file
    writeVideo(vidObj, im2frame(print(fig,'-RGBImage')));
end


% Close video file
close(vidObj);
toc