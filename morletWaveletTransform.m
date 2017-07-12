function [cfx,cfreqso] = morletWaveletTransform(x, fs, cfreqs, morletParam, dim, plotFlag)
% Takes the complex Morlet wavelet transform of time series data and plots
% spectrogram
% 
% ARGUMENTS:
%       x -- 1xT vector time series or multi-dimensional vector with time
%           specified in dimension DIM
%       fs -- Sampling frequency (in Hz)
%       cfreqs -- Vector of centre frequencies
%       morletParam -- (Optional) Morlet wavelet parameter, allows trade
%           between time and frequency resolution (higher is better
%           frequency resolution). Default value 7
%       dim -- (Optional) Specify the dimension of time in X
%       plotFlag -- (Optional) morletWaveletTransform(..., 'plot') means
%          function will create time-frequency plots, with or without
%          morletParam specified. Plots will be averaged across channels.
%          Default is no plotting
%
% OUTPUTS:
%       cfs -- FxT matrix of complex Morlet wavelet coefficients, where F
%           is the number of centre frequencies. If X is multi-dimensional,
%           CFS will be Fx(size(X))
%       cfreqso -- Vector of centre frequencies. Should always be equal to
%           input argument cfreqs, but this can be used to check that CWTFT
%           is not changing the input frequencies.
%
% USAGE:
%{
        fs = 1000;
        t = 0:(1/Fs):2;
        x = chirp(t,3,1,8,'quadratic');
        cfreqs = linspace(1, 10, 100);
        % Use built-in spectrogram plotting with default Morlet parameter
        morletWaveletTransform(x, fs, cfreqs, 'plot');
        % Manually plot spectrogram with smaller parameter
        [cfx, cfreqs] = morletWaveletTransform(x, fs, cfreqs, 3);
        figure
        imagesc(t, cfreqs, abs(cfx))
        axis xy
%}
%
% Rory Townsend, Oct 2016
% rory.townsend@sydney.edu.au


% Set up paramters for morlet wavelet transform
if exist('morletParam', 'var') && strcmp(morletParam, 'plot')
    plotFlag = 'plot';
end
if ~exist('morletParam', 'var') || strcmp(morletParam, 'plot')
    morletParam = 7;
end
if ~isvector(x) && (~exist('dim', 'var') || strcmp(morletParam, 'plot'))
    dim = 2;
end

if ~isvector(x)
    % Reshape input so that time in in the second dimension and other 
    % dimensions are combined
    permOrder = [dim, 1:dim-1, dim+1:ndims(x)];
    x = permute(x, permOrder);
    sx = size(x);
    x = x(:,:);
end

dt = 1/fs;
morletFourierFactor = 4*pi/(morletParam+sqrt(2+morletParam^2));

% Set up a larger output matrix if X has multiple channels / trials
if ~isvector(x)
    cfx = zeros([length(cfreqs), size(x)]);
end

% Old code using CWT rather than CWTFT
% wname = 'morl';
% scales = centfrq(wname)./(cfreqs*dt);

% Set up structure defining scales between min and max pseudo-frequencies
scales = 1./(morletFourierFactor * cfreqs);

% Calculate wavelet coefficients for each channel
if ~isvector(x)
    for ichan = 1:size(x,2)
        %icfs = cwt(x(ichan,:), scales, wname);
        cfstruct = cwtft({x(:,ichan),dt},'scales',scales,'wavelet',{'morl', morletParam});
        cfx(:, :, ichan) = cfstruct.cfs;
    end
    plotVal = squeeze(mean(cfx, 2));
else
    %cfs = cwt(x, scales, wname);
    cfstruct = cwtft({x,dt},'scales',scales,'wavelet',{'morl',morletParam});
    cfx = cfstruct.cfs;
    plotVal = cfx;
end

sc = cfstruct.scales;
cfreqso = 1./(sc*morletFourierFactor);

% Plot data 
if exist('plotFlag', 'var') && (strcmp(plotFlag, 'plot'))
    figure
    % Generate time-frequency power plot
    time = (1:length(x))/fs;
    imagesc(time, cfreqso, abs(plotVal))
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('Morlet wavelet power')
    axis xy
    
    % Generate time-frequency phase plot
    figure
    imagesc(time, cfreqso, angle(plotVal))
    % Note: hsv is MATLAB's only built-in colormap suitable for circular
    % data, but it is very ugly
    if exist('pmkmp_new', 'file') == 2
        colormap(pmkmp_new(256, 'ostwald_o'));
    else
        colormap(hsv)
    end
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('Morlet wavelet phase')
    axis xy
end

% Reshape output to input size
if ~isvector(x)
    cfx = reshape(cfx, [length(cfreqs), sx]);
    cfx = ipermute(cfx, [1, 1+permOrder]);
end


