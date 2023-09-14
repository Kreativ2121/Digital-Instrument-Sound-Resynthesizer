function [fletcherAndMunson40dBEqualLoudnessCurve] = fletcher_and_munson_40dB_curve_generator(frequency)
%FLETCHER_AND_MUNSON_40DB_CURVE_GENERATOR Function that generates a 40dB
%equal loudness curve.
%   This function makes an equal loudness curve, similar to that defined by
%   Fletcher & Munson. Similar, because it is just an approximation, to
%   save computing power.
    fletcherAndMunson40dBEqualLoudnessCurveX = double(.05 + ...
        double(4000./frequency));
    fletcherAndMunson40dBEqualLoudnessCurve = ...
        fletcherAndMunson40dBEqualLoudnessCurveX .* ...
        10.^(-fletcherAndMunson40dBEqualLoudnessCurveX);
end

