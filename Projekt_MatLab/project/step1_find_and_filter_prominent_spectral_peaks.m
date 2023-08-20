function [frequencyPeaks, frequencyPeaksdBFiltered] = ...
    step1_find_and_filter_prominent_spectral_peaks(magnitudeDecibels, ...
    frequency, minimumPeakHeightApprox, minimumPeakHeightLocal, ...
    minimumPeakHeightGlobal, amplitudeRangeLow, frequencyRangeLow, ...
    frequencyRangeHigh)
%step1_find_prominent_spectral_peaks A part of the algorythm that is
%responsible for finding peaks.
%   This part of code finds all peaks in a post STFT data. Then filters all
%   minor or non-important peaks in accordance to set conditions.
    counter = 1;
    frequencyPeaks = [];
    
    % Find peaks
    % Loop on time frames
    for j=1:size(magnitudeDecibels,2)
        % Loop on frequency bins
        for i=2:size(magnitudeDecibels,1)-1
            if(magnitudeDecibels(i,j) >= magnitudeDecibels(i-1,j) && ...
                    magnitudeDecibels(i,j) >= magnitudeDecibels(i+1,j))
                frequencyPeaks(1,counter) = magnitudeDecibels(i,j);
                frequencyPeaks(2,counter) = i; % No. of bin
                frequencyPeaks(3,counter) = j; % No. of time frame
                counter = counter + 1;
            end
        end
    end
    
    % Find valleys for all peaks
    peakValleys = [];
    for i=1:size(frequencyPeaks,2)
        prevValue = frequencyPeaks(1,i);
        for j = frequencyPeaks(2,i):-1:1
            if(magnitudeDecibels(j,frequencyPeaks(3,i)) <= prevValue)
                prevValue = magnitudeDecibels(j,frequencyPeaks(3,i));
            else
                peakValleys(i,1) = prevValue;
                break;
            end
        end
        peakValleys(i,1) = prevValue;
    
        prevValue = frequencyPeaks(1,i);
        for j = frequencyPeaks(2,i):size(magnitudeDecibels,1)
            if(magnitudeDecibels(j,frequencyPeaks(3,i)) <= prevValue)
                prevValue = magnitudeDecibels(j,frequencyPeaks(3,i));
            else
                peakValleys(i,2) = prevValue;
                break;
            end
        end
        peakValleys(i,2) = prevValue;
    end
    peakValleys = peakValleys';
    
    % Assigning Peak Heights and unique index - usable for peak interpolation
    for i=1:size(frequencyPeaks,2)
        frequencyPeaks(4,i) = frequencyPeaks(1,i) - (peakValleys(1,i) + ...
            peakValleys(2,i))/2; % Peak heights in dB
        frequencyPeaks(6,i) = i;
    end
    
    % % Adding magnitude data relative to maxdB (the "negative dB scale")
    maxdB = max(frequencyPeaks(1,:));
    
    counter = 1;
    for i=1:size(frequencyPeaks,2)
        frequencyPeaks(5,counter) = frequencyPeaks(1,i)-maxdB;
        counter = counter+1;
    end
    
    % Find maxAmplitudes in all Time Frames -> useful in discarding small peaks
    MaxAmplitudesInTimeFrames = zeros(1,size(magnitudeDecibels,2));
    for i=1:size(magnitudeDecibels,2)
        MaxAmplitudesInTimeFrames(i) = max(magnitudeDecibels(:,i)) - maxdB;
    end
    
    % Filtering Small Peaks
    frequencyPeaksHeightFiltered = [];
    counter = 1;
    for i=1:size(frequencyPeaks,2)
        maxAmplitudeInTimeFrame = MaxAmplitudesInTimeFrames(frequencyPeaks(3,i));
        minimumPeakHeightLocalThreshold = maxAmplitudeInTimeFrame ...
        +minimumPeakHeightLocal;
    
        % Filtrowanie zarówno na bazie obszaru w którym się znajduje amplituda,
        % jak i na bazie samej wielkości amplitudy pomniejszonej o doliny.
        if((frequencyPeaks(5,i) >= minimumPeakHeightGlobal || ...
                frequencyPeaks(5,i) >= minimumPeakHeightLocalThreshold) ...
                && (frequencyPeaks(4,i) > minimumPeakHeightApprox))
            frequencyPeaksHeightFiltered(1,counter) = frequencyPeaks(1,i);
            frequencyPeaksHeightFiltered(2,counter) = frequencyPeaks(2,i);
            frequencyPeaksHeightFiltered(3,counter) = frequencyPeaks(3,i);
            frequencyPeaksHeightFiltered(4,counter) = frequencyPeaks(4,i);
            frequencyPeaksHeightFiltered(5,counter) = frequencyPeaks(5,i);
            frequencyPeaksHeightFiltered(6,counter) = frequencyPeaks(6,i);
            counter = counter+1;
        end
    end
    
    % Filtering peaks out of the audible frequency range
    frequencyPeaksRangeFiltered = [];
    counter = 1;
    for i=1:size(frequencyPeaksHeightFiltered,2)
        if(frequency(frequencyPeaksHeightFiltered(2,i)) >= frequencyRangeLow ...
                && frequency(frequencyPeaksHeightFiltered(2,i)) <= frequencyRangeHigh)
            frequencyPeaksRangeFiltered(1,counter) = frequencyPeaksHeightFiltered(1,i);
            frequencyPeaksRangeFiltered(2,counter) = frequencyPeaksHeightFiltered(2,i);
            frequencyPeaksRangeFiltered(3,counter) = frequencyPeaksHeightFiltered(3,i);
            frequencyPeaksRangeFiltered(4,counter) = frequencyPeaksHeightFiltered(4,i);
            frequencyPeaksRangeFiltered(5,counter) = frequencyPeaksHeightFiltered(5,i);
            frequencyPeaksRangeFiltered(6,counter) = frequencyPeaksHeightFiltered(6,i);
            counter = counter+1;
        end
    end
    
    % Discard peaks with very low magnitude in general-dB-range
    frequencyPeaksdBFiltered = [];
    counter = 1;
    for i=1:size(frequencyPeaksRangeFiltered,2)
        if(frequencyPeaksRangeFiltered(4,i) >= amplitudeRangeLow)
            frequencyPeaksdBFiltered(1,counter) = frequencyPeaksRangeFiltered(1,i);
            frequencyPeaksdBFiltered(2,counter) = frequencyPeaksRangeFiltered(2,i);
            frequencyPeaksdBFiltered(3,counter) = frequencyPeaksRangeFiltered(3,i);
            frequencyPeaksdBFiltered(4,counter) = frequencyPeaksRangeFiltered(4,i);
            frequencyPeaksdBFiltered(5,counter) = frequencyPeaksRangeFiltered(5,i);
            frequencyPeaksdBFiltered(6,counter) = frequencyPeaksRangeFiltered(6,i);
            counter = counter+1;
        end
    end
end

