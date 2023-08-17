function [phase_inst] = phase_calculation_interpolation(phase_last, phase_now, freq_last, freq_now, m, S)
%PHASE_CALCULATION_ Function calculating a cubic polynomial for phase
%interpolation during resynthesis.
%   A function that will calculate a function given by the formula: θ(m) =
%   ζ + κm + ηm2 + ιm3, which equals to ˆϕ(l−1) + ˆω(l−1)m + ηm2 + ιm3

% 3.25
x = (1/(2*pi)) * ((phase_last + freq_last*S - phase_now) + (S/2)*(freq_now-freq_last));
M = round(x);

% 3.24
eta = (3/power(S,2))*(phase_now - phase_last - freq_last*S + 2*pi*M)-(1/S)*(freq_now-freq_last);
jota = (-2/power(S,3))*(phase_now - phase_last - freq_last*S + 2*pi*M)+(1/power(S,2))*(freq_now-freq_last);

% 3.23
phase_inst = phase_last + freq_last*m + eta*power(m,2) + jota*power(m,3);
% phase_inst = wrapToPi(phase_inst);
end

