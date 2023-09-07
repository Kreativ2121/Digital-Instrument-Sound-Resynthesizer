function [peaks] = normalize_amplitudes(peaks, windowLength, beta)
% normalize_amplitudes Function that normalizes the peak amplitudes in
% regard to the windowLength. This version works only for Keiser window.
    % Normalization scale factor calculation for Kaiser window
    % W0_B = sum(kaiser(windowLength,beta));
    % alpha = double(((windowLength)/128+1)/W0_B);
    ks = kaiser(windowLength,beta);
    ks = ks(floor(length(ks)/2):end);
    W0_B = sum(ks);
    % alpha = double(2/W0_B);
    alpha = double(2/W0_B);
    
    % Normalized amplitude              
    peaks(11,:) = peaks(11,:) * alpha;
end

