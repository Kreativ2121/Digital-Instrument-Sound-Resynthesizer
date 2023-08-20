clc;
clear;
close all;

%% SETTINGS
fftLength = 256;
windowLength = fftLength;
overlapLength = floor(0.75 * windowLength);
hopsize = fftLength - overlapLength;
beta = 6.0;

minimumPeakHeightApprox = 10; %in dB
minimumPeakHeightGlobal = -40; %in dB
minimumPeakHeightLocal = -30; %in dB
frequencyRangeLow = 20; %in hz
frequencyRangeHigh = 16000; %in hz
amplitudeRangeLow = -70; %in dB
amplitudeRangeHigh = 0; %in dB

% MaximumPeakDeviation = 500;
MaximumPeakDeviation = 30; %Większa granica -> mniej trajektorii

%% Read Wave File
% filetitle = "../src/generated/mono/square2000.wav";
% filetitle = "../src/generated/mono/square440.wav";
% filetitle = "../src/generated/mono/square689.wav";
% filetitle = "../src/generated/mono/square2411.wav";
filetitle = "../src/generated/mono/sine440.wav";
% filetitle = "../src/generated/mono/sine689.wav";
% filetitle = "../src/generated/mono/saw689.wav";
% filetitle = "../src/generated/mono/chirp440_2000.wav";
% filetitle = "../src/generated/mono/chirp2000_8000.wav";
% filetitle = "../src/generated/mono/sine2000.wav";
% filetitle = "../src/generated/mono/square2000_additivesynthesis.wav";
% filetitle = "../src/download/CantinaBand3.wav";

%% STEP 1
[audioData, fs] = audioread(filetitle);
auInfo = audioinfo(filetitle);

[magnitude,frequency,time] = stft(audioData,fs, ...
    Window=kaiser(windowLength,beta), FFTLength=fftLength, ...
    OverlapLength=overlapLength, FrequencyRange="onesided");
% plot_wavelet_and_stft(audioData, fs, fftLength, windowLength, ...
%     overlapLength, beta);

magnitudeDecibels = 20*log10(abs(magnitude));

[frequencyPeaks, frequencyPeaksdBFiltered] = step1_find_and_filter_prominent_spectral_peaks( ...
    magnitudeDecibels, frequency, minimumPeakHeightApprox, ...
    minimumPeakHeightLocal, minimumPeakHeightGlobal, amplitudeRangeLow, ...
    frequencyRangeLow, frequencyRangeHigh);

% Creating equal loudness curve
fletcher_and_mundson_40dB = fletcher_and_mundson_40dB_curve_generator(frequency);

%% STEP 2
peaks = step2_peak_detection_and_interpolation(frequencyPeaksdBFiltered, ...
    magnitude, magnitudeDecibels, frequency, fletcher_and_mundson_40dB);
peaks = normalize_amplitudes(peaks, windowLength, beta);

%% STEP 3
Trajectories = step3_assign_peak_frequency_trajectories(peaks, MaximumPeakDeviation);

%% STEP 4
audioDataLength = length(audioData);
output = step4_resynthesize(Trajectories, fs, hopsize, audioDataLength);

audiowrite("output.wav",output,fs);

%% COMPARISONS
plot_results(audioData, output, fs, windowLength, beta, fftLength, overlapLength)










