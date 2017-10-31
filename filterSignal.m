function filteredSignal = filterSignal(x, fLow, fHigh, Fs)
% FILTERSIGNAL band-pass filters the signal X with the desired cutoff
% frequencies FLOW and FHIGH, given the sampling frequency FS.
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

N = 8;  % Filter order

if nargin==4 && fLow>0 && fHigh>fLow
    % Design Butterworth band-pass filter
    h = fdesign.bandpass('N,F3dB1,F3dB2',N,fLow,fHigh,Fs);
    Hd = design(h, 'butter');
    set(Hd, 'Arithmetic', 'double');
    
    SOS = Hd.sosMatrix;
    G = Hd.ScaleValues;
    % Filter signals forwards and backwards to avoid phase distortion
    filteredSignal = filtfilt(SOS,G,x);

else
    % Throw an error if inputs are invalid
    error('Invalid inputs for filtering!');
end

end