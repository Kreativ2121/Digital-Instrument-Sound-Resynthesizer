clc;
clear;
xLimitation = [duration(0,0,0,0) duration(0,0,0,100)];

%% Read Wave File
filetitle = "src/generated/mono/square2000.wav";
[audioData,fs] = audioread(filetitle);
auInfo = audioinfo(filetitle);

% Zapisanie do zmiennej tylko nazwy pliku ze ścieżki
file = strsplit(filetitle,'/');
file = file(end);

%% Plot wavelet
t = seconds(0:1/fs:(size(audioData,1)-1)/fs);

subplot(2,1,1)
plot(t,audioData)
title("Przebieg")
xlabel("Czas")
ylabel("Amplituda")
legend("Kanał 1", "Kanał 2")
xlim(xLimitation)
%xlim("tight")
ylim([-1 1])

%% ALGORITHM
subplot(2,1,2)

%%STEP 1
stft(audioData,fs)
[magnitude,frequency,time] = stft(audioData,fs);
% MagnitudeDecibels = mag2db(abs(magnitude));

%%STEP 2
Ray = abs(magnitude);
PointAmplitude = deg2rad(magnitude);

%%STEP 3
MagnitudeDecibels = 20*log10(Ray);

%%STEP 4
MinimumPeakHeight = 15; %in dB
FrequencyRangeLow = 20; %in hz
FrequencyRangeHigh = 16000; %in hz
AmplitudeRangeLow = -70; %in dB ?
AmplitudeRangeHigh = 0; %in dB ?

counter = 1;
FrequencyPeaks = [];

% Find peaks
% Loop on time frames
% TODO Sprawić by krawędziowe elementy mogły być wliczane jako piki
for j=1:size(MagnitudeDecibels,2)
    % Loop on frequency bins
    for i=2:size(MagnitudeDecibels,1)-1
        if(MagnitudeDecibels(i,j) >= MagnitudeDecibels(i-1,j) && MagnitudeDecibels(i,j) >= MagnitudeDecibels(i+1,j))
            FrequencyPeaks(1,counter) = MagnitudeDecibels(i,j);
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
        if(MagnitudeDecibels(j,FrequencyPeaks(3,i)) <= prevValue)
            prevValue = MagnitudeDecibels(j,FrequencyPeaks(3,i));
        else
            PeakValleys(i,1) = prevValue;
            break;
        end
    end

    prevValue = FrequencyPeaks(1,i);
    for j = FrequencyPeaks(2,i):size(MagnitudeDecibels,1)
        if(MagnitudeDecibels(j,FrequencyPeaks(3,i)) <= prevValue)
            prevValue = MagnitudeDecibels(j,FrequencyPeaks(3,i));
        else
            PeakValleys(i,2) = prevValue;
            break;
        end
    end
end
PeakValleys = PeakValleys';

% Assigning Peak Heights
for i=1:size(FrequencyPeaks,2)
    FrequencyPeaks(4,i) = FrequencyPeaks(1,i) - (PeakValleys(1,i) + PeakValleys(2,i))/2; % Peak heights in dB
end

% Filtering Small Peaks
FrequencyPeaksHeightFiltered = [];
counter = 1;
for i=1:size(FrequencyPeaks,2)
    if(FrequencyPeaks(4,i) >= MinimumPeakHeight)
        FrequencyPeaksHeightFiltered(1,counter) = FrequencyPeaks(1,i);
        FrequencyPeaksHeightFiltered(2,counter) = FrequencyPeaks(2,i);
        FrequencyPeaksHeightFiltered(3,counter) = FrequencyPeaks(3,i);
        FrequencyPeaksHeightFiltered(4,counter) = FrequencyPeaks(4,i);
        counter = counter+1;
    end
end

% Obtain frequency bins
kl = FrequencyRangeLow*size(frequency,1)/fs;
kh = FrequencyRangeHigh*size(frequency,1)/fs;

% Filtering peaks out of the audible frequency range
FrequencyPeaksRangeFiltered = [];
counter = 1;
for i=1:size(FrequencyPeaksHeightFiltered,2)
    if(frequency(FrequencyPeaksHeightFiltered(2,i)) >= FrequencyRangeLow && frequency(FrequencyPeaksHeightFiltered(2,i)) <= FrequencyRangeHigh)
        FrequencyPeaksRangeFiltered(1,counter) = FrequencyPeaksHeightFiltered(1,i);
        FrequencyPeaksRangeFiltered(2,counter) = FrequencyPeaksHeightFiltered(2,i);
        FrequencyPeaksRangeFiltered(3,counter) = FrequencyPeaksHeightFiltered(3,i);
        FrequencyPeaksRangeFiltered(4,counter) = FrequencyPeaksHeightFiltered(4,i);
        counter = counter+1;
    end
end

% Adding a row that will show values relative to max dB in whole sound
maxdB = max(FrequencyPeaksRangeFiltered(4,:));

counter = 1;
for i=1:size(FrequencyPeaksRangeFiltered,2)
    FrequencyPeaksRangeFiltered(5,counter) = FrequencyPeaksHeightFiltered(4,i)-maxdB;
    counter = counter+1;
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
        counter = counter+1;
    end
end

% Creating equal loudness curve
FletcherMundson40LikeLoudnessCurveX = double(.05 + double(4000./frequency));
FletcherMundson40LikeLoudnessCurveX = FletcherMundson40LikeLoudnessCurveX(ceil(length(FletcherMundson40LikeLoudnessCurveX)/2):end);
FletcherMundson40LikeLoudnessCurve = FletcherMundson40LikeLoudnessCurveX .* 10.^(-FletcherMundson40LikeLoudnessCurveX);
% plot(frequency(ceil(length(frequency)/2):end),FletcherMundson40LikeLoudnessCurve)

% Choosing only positive frequency values - necessary?
FrequencyPeaksPositive = [];
counter = 1;
for i=1:size(FrequencyPeaksdBFiltered,2)
    if(FrequencyPeaksdBFiltered(2,i) >= ceil(length(frequency)/2))
        FrequencyPeaksPositive(1,counter) = FrequencyPeaksdBFiltered(1,i);
        FrequencyPeaksPositive(2,counter) = FrequencyPeaksdBFiltered(2,i);
        FrequencyPeaksPositive(3,counter) = FrequencyPeaksdBFiltered(3,i);
        FrequencyPeaksPositive(4,counter) = FrequencyPeaksdBFiltered(4,i);
        FrequencyPeaksPositive(5,counter) = FrequencyPeaksdBFiltered(5,i);
        counter = counter+1;
    end
end

% Applying equal loudness curve on magnitude as a sixth row
% for i=1:size(FrequencyPeaksRangeFiltered,2)
%     FrequencyPeaksdBFiltered(5,counter) = FrequencyPeaksRangeFiltered(5,i) - FletcherMundson40LikeLoudnessCurve(i);
%     counter = counter+1;
% end

%%STEP 5

%%STEP 6

%%STEP 7

%%STEP 8 - RESYNTHESIS





%% NOTES

% Bez argumentów - wyświetla wykres
% stft(audioData,fs,"Window",gausswin(size(audioData,1)/500))

%%% Manualne wyświetlenie wykresu 
% sdb = mag2db(abs(s));
% mesh(t,f/1000,sdb);
% 
% cc = max(sdb(:))+[-60 0];
% ax = gca;
% ax.CLim = cc;
% view(2)
% colorbar
% title("Short Time Fourier Transform (STFT)")
% xlabel("Time (s)")
% ylabel("Frequency (kHz)")
