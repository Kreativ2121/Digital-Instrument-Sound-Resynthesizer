function [peaks] = normalize_amplitudes(peaks, windowLength, beta)
% normalize_amplitudes Function that normalizes the peak amplitudes in
% regard to the windowLength. This version works only for Keiser window.
    % Normalization scale factor calculation for Kaiser window
    W0_B = sum(kaiser(196,beta));
    alpha = double(2/W0_B);
    
    % Normalized amplitude
    N_Amp = peaks(11,:) * alpha;                    
    peaks(11,:) = N_Amp;
end

