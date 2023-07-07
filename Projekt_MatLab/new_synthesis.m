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
MinimumPeakHeight = 3; %in dB
FrequencyRangeLow = 20; %in hz
FrequencyRangeHigh = 16000; %in hz
AmplitudeRangeLow = 0; %in dB ?
AmplitudeRangeHigh = 80; %in dB ?

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
