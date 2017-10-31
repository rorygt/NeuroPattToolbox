function fx = filterSignal(x, fLow, fHigh, Fs)
% FILTERSIGNAL band-pass filters the signal X with the desired cutoff
% frequencies FLOW and FHIGH, given the sampling frequency FS.
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

N = 8;  % Filter order

if ~isvector(x)
    % Reshape input so that time is in the first dimension and other 
    % dimensions are combined
    permOrder = [dim, 1:dim-1, dim+1:ndims(x)];
    x = permute(x, permOrder);
    sx = size(x);
    x = x(:,:);
end

fx = zeros(size(x));

if nargin==4 && fLow>0 && fHigh>fLow
    % Design Butterworth band-pass filter
    h = fdesign.bandpass('N,F3dB1,F3dB2',N,fLow,fHigh,Fs);
    Hd = design(h, 'butter');
    set(Hd, 'Arithmetic', 'double');
    
    SOS = Hd.sosMatrix;
    G = Hd.ScaleValues;
    
    for ii = 1:size(x,2)
        % Filter signals forwards and backwards to avoid phase distortion
        ifx = filtfilt(SOS,G,x(:,ii));
        fx(:,ii) = ifx;
    end
else
    % Throw an error if inputs are invalid
    error('Invalid inputs for filtering!');
end

% Reshape output to input size
if ~isvector(x)
    fx = ipermute(reshape(fx, sx), permorder);
end

end