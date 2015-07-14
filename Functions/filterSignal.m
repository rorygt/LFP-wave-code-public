function filteredSignal = filterSignal(LFP, fLow, fHigh, Fs)

% Filters the signal (LFP) with the desired cutoff frequencies (fLow and
% fHigh), given the sampling frequency (Fs, should be 1024 for Utah data).

N = 8;  % Filter order

if fLow ~= 0
    
    % Design band pass filter
    h = fdesign.bandpass('N,F3dB1,F3dB2',N,fLow,fHigh,Fs);
    Hd = design(h, 'butter');
    set(Hd, 'Arithmetic', 'double');
    
    % Filter forwards and backwards in time using the filter defined by
    % second order section matrix SOS and scale values G
    SOS = Hd.sosMatrix;
    G = Hd.ScaleValues;
    filteredSignal = filtfilt(SOS,G,LFP);

else
    % If fLow == 0, the fdesign function will not work, so just return the
    % unfiltered data,
    display('Lower cutoff frequency must not be zero!');
    filteredSignal = LFP;
    return
    
end

end