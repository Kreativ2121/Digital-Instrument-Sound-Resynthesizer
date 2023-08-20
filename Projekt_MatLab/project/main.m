clc;
clear;
close all;

%% SETTINGS
fftlength = 256;
windowlength = fftlength;
overlaplength = floor(0.75 * windowlength);
hopsize = fftlength - overlaplength;
beta = 6.0;

MinimumPeakHeightApprox = 10; %in dB
MinimumPeakHeightGlobal = -40; %in dB
MinimumPeakHeightLocal = -30; %in dB
FrequencyRangeLow = 20; %in hz
FrequencyRangeHigh = 16000; %in hz
AmplitudeRangeLow = -70; %in dB
AmplitudeRangeHigh = 0; %in dB

% MaximumPeakDeviation = 500;
MaximumPeakDeviation = 30; %Większa granica -> mniej trajektorii

%% Read Wave File
% filetitle = "src/generated/mono/square2000.wav";
% filetitle = "src/generated/mono/square440.wav";
% filetitle = "src/generated/mono/square689.wav";
% filetitle = "src/generated/mono/square2411.wav";
filetitle = "../src/generated/mono/sine440.wav";
% filetitle = "src/generated/mono/sine689.wav";
%filetitle = "src/generated/mono/saw689.wav";
% filetitle = "src/generated/mono/chirp440_2000.wav";
% filetitle = "src/generated/mono/chirp2000_8000.wav";
% filetitle = "src/generated/mono/sine2000.wav";
% filetitle = "src/generated/mono/square2000_additivesynthesis.wav";
% filetitle = "src/download/CantinaBand3.wav";

%% STEP 1
[audioData,fs] = audioread(filetitle);
auInfo = audioinfo(filetitle);

[magnitude,frequency,time] = stft(audioData,fs, Window=kaiser(windowlength,beta), FFTLength=fftlength, OverlapLength=overlaplength, FrequencyRange="onesided");
% plot_wavelet_and_stft(audioData, fs, fftlength, windowlength, overlaplength, beta);

magnitudeDecibels = 20*log10(abs(magnitude));

[FrequencyPeaks, FrequencyPeaksdBFiltered] = step1_find_and_filter_prominent_spectral_peaks(magnitudeDecibels, frequency, MinimumPeakHeightApprox, MinimumPeakHeightLocal, MinimumPeakHeightGlobal, AmplitudeRangeLow, FrequencyRangeLow, FrequencyRangeHigh);

% Creating equal loudness curve
fletcher_and_mundson_40dB = fletcher_and_mundson_40dB_curve_generator(frequency);

%% STEP 2
Peaks = step2_peak_detection_and_interpolation(FrequencyPeaksdBFiltered, magnitude, magnitudeDecibels, frequency, fletcher_and_mundson_40dB);
Peaks = normalize_amplitudes(Peaks, windowlength, beta);

%% STEP 3
Trajectories = step3_assign_peak_frequency_trajectories(Peaks, MaximumPeakDeviation);

%% STEP 4
audioDataLength = length(audioData);
output = step4_resynthesize(Trajectories, fs, hopsize, audioDataLength);

audiowrite("output.wav",output,fs);

%% COMPARISONS
% plot_results(audioData, output, fs, windowlength, beta, fftlength, overlaplength)










