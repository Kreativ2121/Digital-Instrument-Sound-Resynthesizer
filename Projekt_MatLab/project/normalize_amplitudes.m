function [Peaks] = normalize_amplitudes(Peaks, windowlength, beta)
% normalize_amplitudes Summary of this function goes here
%   Detailed explanation goes here
    % Normalization scale factor calculation for Kaiser window
    W0_B = sum(kaiser(windowlength,beta));
    alpha = double(2/W0_B);
    
    % Normalized amplitude
    N_Amp = Peaks(11,:) * alpha;                    
    Peaks(11,:) = N_Amp;
end

