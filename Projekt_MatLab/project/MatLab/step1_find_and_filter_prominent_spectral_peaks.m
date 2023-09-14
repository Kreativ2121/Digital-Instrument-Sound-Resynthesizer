function [frequencyPeaks, frequencyPeaksFiltered] = ...
    step1_find_and_filter_prominent_spectral_peaks(magnitudeDecibels, ...
    frequency, minimumPeakHeightLocal, ...
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
                    magnitudeDecibels(i,j) >= magnitudeDecibels(i+1,j) && ...
                    ~isinf(magnitudeDecibels(i,j)))
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
    
    % Adding magnitude data relative to maxdB (the "negative dB scale")
    maxdB = max(frequencyPeaks(1,:));

    % Find maxAmplitudes in all Time Frames -> useful in discarding small peaks
    MaxAmplitudesInTimeFrames = zeros(1,size(magnitudeDecibels,2));
    for i=1:size(magnitudeDecibels,2)
        MaxAmplitudesInTimeFrames(i) = max(magnitudeDecibels(:,i)) - maxdB;
    end

    % Assigning Peak Heights and unique index - usable for peak interpolation
    for i=1:size(frequencyPeaks,2)
        frequencyPeaks(4,i) = frequencyPeaks(1,i) - (peakValleys(1,i) + ...
            peakValleys(2,i))/2; % Peak heights in dB
        frequencyPeaks(5,i) = frequencyPeaks(1,i)-maxdB;
        frequencyPeaks(6,i) = i;
    end
    
    % Filtering Small Peaks
    frequencyPeaksFiltered = [];
    counter = 1;
    for i=1:size(frequencyPeaks,2)
        maxAmplitudeInTimeFrame = MaxAmplitudesInTimeFrames(frequencyPeaks(3,i));
        minimumPeakHeightLocalThreshold = maxAmplitudeInTimeFrame ...
        +minimumPeakHeightLocal;
    
        % Filtering based on an amplitude height
        % Filtering peaks out of the audible frequency range
        % Discard peaks with very low magnitude in general-dB-range
        if((frequencyPeaks(5,i) >= minimumPeakHeightGlobal ...
                || frequencyPeaks(5,i) >= minimumPeakHeightLocalThreshold) ...
                && frequency(frequencyPeaks(2,i)) >= frequencyRangeLow ...
                && frequency(frequencyPeaks(2,i)) <= frequencyRangeHigh ...
                && frequencyPeaks(4,i) >= amplitudeRangeLow)
            frequencyPeaksFiltered(1,counter) = frequencyPeaks(1,i);
            frequencyPeaksFiltered(2,counter) = frequencyPeaks(2,i);
            frequencyPeaksFiltered(3,counter) = frequencyPeaks(3,i);
            frequencyPeaksFiltered(4,counter) = frequencyPeaks(4,i);
            frequencyPeaksFiltered(5,counter) = frequencyPeaks(5,i);
            frequencyPeaksFiltered(6,counter) = frequencyPeaks(6,i);
            counter = counter+1;
        end
    end
end

