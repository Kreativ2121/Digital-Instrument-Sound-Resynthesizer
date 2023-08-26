function [peaks] = step2_peak_detection_and_interpolation(frequencyPeaksdBFiltered, magnitudeDecibels, frequency, fletcher_and_mundson_40dB)
%step2_peak_detection_and_interpolation A part of the algorithm that is
%responsible for peak interpolation.
%   this part of the algorythm is responsible for peak interpolation. Based
%   on a previously calculated data it tries to correct the location of the
%   peak, because the data after STFT is not precise enough.
    peaks = [];
    distanceBetweenFreqBins = frequency(2)-frequency(1);
    counter = 1;
    for i=1:size(frequencyPeaksdBFiltered,2)
        peaks(1,counter) = frequencyPeaksdBFiltered(1,i);
        peaks(2,counter) = frequencyPeaksdBFiltered(2,i);
        peaks(3,counter) = frequencyPeaksdBFiltered(3,i);
        peaks(4,counter) = frequencyPeaksdBFiltered(4,i);
        peaks(5,counter) = frequencyPeaksdBFiltered(5,i);
        peaks(6,counter) = frequencyPeaksdBFiltered(6,i);
    
        %%ALPHA
        if(frequencyPeaksdBFiltered(2,i)-1 == 0)
            peaks(7,counter) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i), ...
                frequencyPeaksdBFiltered(3,i)) ...
                -fletcher_and_mundson_40dB(frequencyPeaksdBFiltered(2,i)));
        else
            % Zakładamy, że valley jest takie samo dla elementów znajdujących się obok pików
            peaks(7,counter) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i)-1, ...
                frequencyPeaksdBFiltered(3,i)) ...
                -fletcher_and_mundson_40dB(frequencyPeaksdBFiltered(2,i)));
        end
    
        %%BETA
        peaks(8,counter) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i), ...
            frequencyPeaksdBFiltered(3,i)) ...
            -fletcher_and_mundson_40dB(frequencyPeaksdBFiltered(2,i)));
    
        %%GAMMA
        if(frequencyPeaksdBFiltered(2,i)+1 == length(frequencyPeaksdBFiltered))
            peaks(9,counter) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i), ...
                frequencyPeaksdBFiltered(3,i)) ...
                -fletcher_and_mundson_40dB(frequencyPeaksdBFiltered(2,i)));
        else
            peaks(9,counter) = (magnitudeDecibels(frequencyPeaksdBFiltered(2,i)+1, ...
                frequencyPeaksdBFiltered(3,i)) ...
                -fletcher_and_mundson_40dB(frequencyPeaksdBFiltered(2,i)));
        end
    
        %Assign parabola peak location
        peaks(10,counter) = ((peaks(7,counter)- ...
            peaks(9,counter))/(peaks(7,counter)-2*peaks(8,counter) ...
            +peaks(9,counter)))/2;
    
        %Assign magnitude peak height
        peaks(11,counter) = peaks(8,counter) - ((peaks(7,counter)- ...
            peaks(9,counter))/4)*peaks(10,counter);
    
        %Assign True Peak Location (in bins)
        peaks(12,counter) = frequency(peaks(2,counter)) + ...
            peaks(10,counter)*distanceBetweenFreqBins;
    
        counter = counter+1;
    end
end

