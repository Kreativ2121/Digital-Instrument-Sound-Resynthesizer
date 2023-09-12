function [peaks] = step2_interpolation_nofm40db(frequencyPeaksdBFiltered, magnitudeDecibels, frequency)
%step2_interpolation A part of the algorithm that is
%responsible for peak interpolation.
%   this part of the algorythm is responsible for peak interpolation. Based
%   on a previously calculated data it tries to correct the location of the
%   peak, because the data after STFT is not precise enough.
    peaks = [];
    distanceBetweenFreqBins = frequency(2)-frequency(1);
    for i=1:size(frequencyPeaksdBFiltered,2)
        peaks(1,i) = frequencyPeaksdBFiltered(1,i);
        peaks(2,i) = frequencyPeaksdBFiltered(2,i);
        peaks(3,i) = frequencyPeaksdBFiltered(3,i);
        peaks(4,i) = frequencyPeaksdBFiltered(4,i);
        peaks(5,i) = frequencyPeaksdBFiltered(5,i);
        peaks(6,i) = frequencyPeaksdBFiltered(6,i);

        %%ALPHA
        if(frequencyPeaksdBFiltered(2,i)-1 == 0)
            peaks(7,i) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i), ...
                frequencyPeaksdBFiltered(3,i)));
        else
            % Zakładamy, że valley jest takie samo dla elementów znajdujących się obok pików
            peaks(7,i) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i)-1, ...
                frequencyPeaksdBFiltered(3,i)));
        end

        %%BETA
        peaks(8,i) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i), ...
            frequencyPeaksdBFiltered(3,i)));

        %%GAMMA
        if(frequencyPeaksdBFiltered(2,i)+1 == length(frequencyPeaksdBFiltered))
            peaks(9,i) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i), ...
                frequencyPeaksdBFiltered(3,i)));
        else
            peaks(9,i) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i)+1, ...
                frequencyPeaksdBFiltered(3,i)));
        end

        %Assign parabola peak location
        peaks(10,i) = ((peaks(7,i)- ...
            peaks(9,i))/(peaks(7,i)-2*peaks(8,i) ...
            +peaks(9,i)))/2;

        %Assign magnitude peak height
        peaks(11,i) = peaks(8,i) - ((peaks(7,i)- ...
            peaks(9,i))/4)*peaks(10,i);

        %Assign True Peak Location (in bins)
        peaks(12,i) = frequency(peaks(2,i)) + ...
            peaks(10,i)*distanceBetweenFreqBins;
    end
end

