function [peaks] = normalize_amplitudes(peaks, windowLength, beta)
% normalize_amplitudes Function that normalizes the peak amplitudes in
% regard to the windowLength. This version works only for Keiser window.
    % Normalization scale factor calculation for Kaiser window
    W0_B = sum(kaiser(windowLength,beta));
    alpha = double(5.05/W0_B);

    % Normalized amplitude              
    peaks(11,:) = peaks(11,:) * alpha;
end