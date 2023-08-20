function [fletcherAndMundson40dBEqualLoudnessCurve] = fletcher_and_mundson_40dB_curve_generator(frequency)
%FLETCHER_AND_MUNDSON_40DB_CURVE_GENERATOR Function that generates a 40dB
%equal loudness curve.
%   This function makes an equal loudness curve, similar to that defined by
%   Fletcher & Mundson. Similar, because it is just an approximation, to
%   save computing power.
    fletcherAndMundson40dBEqualLoudnessCurveX = double(.05 + double(4000./frequency));
    fletcherAndMundson40dBEqualLoudnessCurve = fletcherAndMundson40dBEqualLoudnessCurveX .* 10.^(-fletcherAndMundson40dBEqualLoudnessCurveX);
end

