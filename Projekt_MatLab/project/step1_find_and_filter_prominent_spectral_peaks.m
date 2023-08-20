function [FrequencyPeaks, FrequencyPeaksdBFiltered] = step1_find_and_filter_prominent_spectral_peaks(magnitudeDecibels, frequency, MinimumPeakHeightApprox, MinimumPeakHeightLocal, MinimumPeakHeightGlobal, AmplitudeRangeLow, FrequencyRangeLow, FrequencyRangeHigh)
%step1_find_prominent_spectral_peaks Step 1 of the algorythm
%   Detailed explanation goes here
    counter = 1;
    FrequencyPeaks = [];
    
    % Find peaks
    % Loop on time frames
    for j=1:size(magnitudeDecibels,2)
        % Loop on frequency bins
        for i=2:size(magnitudeDecibels,1)-1
            if(magnitudeDecibels(i,j) >= magnitudeDecibels(i-1,j) && magnitudeDecibels(i,j) >= magnitudeDecibels(i+1,j))
                FrequencyPeaks(1,counter) = magnitudeDecibels(i,j);
                FrequencyPeaks(2,counter) = i; % No. of bin
                FrequencyPeaks(3,counter) = j; % No. of time frame
                counter = counter + 1;
            end
        end
    end
    
    % Find valleys for all peaks
    PeakValleys = [];
    for i=1:size(FrequencyPeaks,2)
        prevValue = FrequencyPeaks(1,i);
        for j = FrequencyPeaks(2,i):-1:1
            if(magnitudeDecibels(j,FrequencyPeaks(3,i)) <= prevValue)
                prevValue = magnitudeDecibels(j,FrequencyPeaks(3,i));
            else
                PeakValleys(i,1) = prevValue;
                break;
            end
        end
        PeakValleys(i,1) = prevValue;
    
        prevValue = FrequencyPeaks(1,i);
        for j = FrequencyPeaks(2,i):size(magnitudeDecibels,1)
            if(magnitudeDecibels(j,FrequencyPeaks(3,i)) <= prevValue)
                prevValue = magnitudeDecibels(j,FrequencyPeaks(3,i));
            else
                PeakValleys(i,2) = prevValue;
                break;
            end
        end
        PeakValleys(i,2) = prevValue;
    end
    PeakValleys = PeakValleys';
    
    % Assigning Peak Heights and unique index - usable for peak interpolation
    for i=1:size(FrequencyPeaks,2)
        FrequencyPeaks(4,i) = FrequencyPeaks(1,i) - (PeakValleys(1,i) + PeakValleys(2,i))/2; % Peak heights in dB
        FrequencyPeaks(6,i) = i;
    end
    
    % % Adding magnitude data relative to maxdB (the "negative dB scale")
    maxdB = max(FrequencyPeaks(1,:));
    
    counter = 1;
    for i=1:size(FrequencyPeaks,2)
        FrequencyPeaks(5,counter) = FrequencyPeaks(1,i)-maxdB;
        counter = counter+1;
    end
    
    % Find maxAmplitudes in all Time Frames -> useful in discarding small peaks
    MaxAmplitudesInTimeFrames = zeros(1,size(magnitudeDecibels,2));
    for i=1:size(magnitudeDecibels,2)
        MaxAmplitudesInTimeFrames(i) = max(magnitudeDecibels(:,i)) - maxdB;
    end
    
    % Filtering Small Peaks
    FrequencyPeaksHeightFiltered = [];
    counter = 1;
    for i=1:size(FrequencyPeaks,2)
        maxAmplitudeInTimeFrame = MaxAmplitudesInTimeFrames(FrequencyPeaks(3,i));
        MinimumPeakHeightLocalThreshold = maxAmplitudeInTimeFrame + MinimumPeakHeightLocal;
    
        % Filtrowanie zarówno na bazie obszaru w którym się znajduje amplituda,
        % jak i na bazie samej wielkości amplitudy pomniejszonej o doliny.
        if((FrequencyPeaks(5,i) >= MinimumPeakHeightGlobal || FrequencyPeaks(5,i) >= MinimumPeakHeightLocalThreshold) && (FrequencyPeaks(4,i) > MinimumPeakHeightApprox))
            FrequencyPeaksHeightFiltered(1,counter) = FrequencyPeaks(1,i);
            FrequencyPeaksHeightFiltered(2,counter) = FrequencyPeaks(2,i);
            FrequencyPeaksHeightFiltered(3,counter) = FrequencyPeaks(3,i);
            FrequencyPeaksHeightFiltered(4,counter) = FrequencyPeaks(4,i);
            FrequencyPeaksHeightFiltered(5,counter) = FrequencyPeaks(5,i);
            FrequencyPeaksHeightFiltered(6,counter) = FrequencyPeaks(6,i);
            counter = counter+1;
        end
    end
    
    % Filtering peaks out of the audible frequency range
    FrequencyPeaksRangeFiltered = [];
    counter = 1;
    for i=1:size(FrequencyPeaksHeightFiltered,2)
        if(frequency(FrequencyPeaksHeightFiltered(2,i)) >= FrequencyRangeLow && frequency(FrequencyPeaksHeightFiltered(2,i)) <= FrequencyRangeHigh)
            FrequencyPeaksRangeFiltered(1,counter) = FrequencyPeaksHeightFiltered(1,i);
            FrequencyPeaksRangeFiltered(2,counter) = FrequencyPeaksHeightFiltered(2,i);
            FrequencyPeaksRangeFiltered(3,counter) = FrequencyPeaksHeightFiltered(3,i);
            FrequencyPeaksRangeFiltered(4,counter) = FrequencyPeaksHeightFiltered(4,i);
            FrequencyPeaksRangeFiltered(5,counter) = FrequencyPeaksHeightFiltered(5,i);
            FrequencyPeaksRangeFiltered(6,counter) = FrequencyPeaksHeightFiltered(6,i);
            counter = counter+1;
        end
    end
    
    % Discard peaks with very low magnitude in general-dB-range
    FrequencyPeaksdBFiltered = [];
    counter = 1;
    for i=1:size(FrequencyPeaksRangeFiltered,2)
        if(FrequencyPeaksRangeFiltered(4,i) >= AmplitudeRangeLow)
            FrequencyPeaksdBFiltered(1,counter) = FrequencyPeaksRangeFiltered(1,i);
            FrequencyPeaksdBFiltered(2,counter) = FrequencyPeaksRangeFiltered(2,i);
            FrequencyPeaksdBFiltered(3,counter) = FrequencyPeaksRangeFiltered(3,i);
            FrequencyPeaksdBFiltered(4,counter) = FrequencyPeaksRangeFiltered(4,i);
            FrequencyPeaksdBFiltered(5,counter) = FrequencyPeaksRangeFiltered(5,i);
            FrequencyPeaksdBFiltered(6,counter) = FrequencyPeaksRangeFiltered(6,i);
            counter = counter+1;
        end
    end
end

