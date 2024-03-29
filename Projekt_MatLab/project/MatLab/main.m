clc;
clear;
close all;

%% SETTINGS
fftLength = 512;
windowLength = 512;
overlapLength = floor(0.75 * windowLength);
hopsize = fftLength - overlapLength;
beta = 6.0;

minimumPeakHeightGlobal = -30; %in dB
minimumPeakHeightLocal = -15; %in dB
frequencyRangeLow = 20; %in hz
frequencyRangeHigh = 16000; %in hz
amplitudeRangeLow = 10; %in dB 

MaximumPeakDeviation = 30; 

%% Read Wave File
% filetitle = "../../src/generated/mono/square2000.wav";
% filetitle = "../../src/generated/mono/square440.wav";
% filetitle = "../../src/generated/mono/square689.wav";
% filetitle = "../../src/generated/mono/square2411.wav";
% filetitle = "../../src/generated/mono/sine440.wav";
% filetitle = "../../src/generated/mono/sine689.wav";
% filetitle = "../../src/generated/mono/saw689.wav";
% filetitle = "../../src/generated/mono/chirp440_2000.wav";
% filetitle = "../../src/generated/mono/chirp2000_8000.wav";
% filetitle = "../../src/generated/mono/chirp2000_6000.wav";
% filetitle = "../../src/generated/mono/chirp2000_14000.wav";
% filetitle = "../../src/generated/mono/chirp14000_2000.wav";
filetitle = "../../src/generated/mono/sine2000.wav";
% filetitle = "../../src/generated/mono/silence.wav";
% filetitle = "../../src/generated/mono/square2000_additivesynthesis.wav";
% filetitle = "../../src/records/kross/mono/KGP_C.wav";
% filetitle = "../../src/generated/mono/2sines4000_1000.wav";
% filetitle = "../../src/generated/mono/2sines2000_1000_2.wav";
% filetitle = "../../src/generated/mono/2sines4000_1000_2.wav";
% filetitle = "../../src/generated/mono/2sines1000_440.wav";
% filetitle = "../../src/generated/mono/2sines4000_2000.wav";
% filetitle = "../../src/generated/mono/2sines2000_3000.wav";
% filetitle = "../../src/generated/mono/2sines4000_2000_ml.wav";
% filetitle = "../../src/generated/mono/4sines_440_1000_3000_8000.wav";
% filetitle = "../../src/generated/mono/4sines_440_1000_3000_1600.wav";
% filetitle = "../../src/records/minilogue/synt2p.wav";
% filetitle = "../../src/records/minilogue/synt1p.wav";
% filetitle = "../../src/generated/mono/silence_then_sound.wav";
% filetitle = "../../src/download/CantinaBand3.wav";

%% STEP 1
[audioData, fs] = audioread(filetitle);
auInfo = audioinfo(filetitle);

[magnitude,frequency,time] = stft(audioData,fs, ...
    Window=kaiser(windowLength,beta), FFTLength=fftLength, ...
    OverlapLength=overlapLength, FrequencyRange="onesided");
plot_waveform_and_stft(audioData, fs, fftLength, windowLength, ...
    overlapLength, beta);

magnitudeDecibels = 20*log10(abs(magnitude));

[frequencyPeaks, frequencyPeaksFiltered] = step1_find_and_filter_prominent_spectral_peaks( ...
    magnitudeDecibels, frequency, minimumPeakHeightLocal, ...
    minimumPeakHeightGlobal, amplitudeRangeLow, ...
    frequencyRangeLow, frequencyRangeHigh);

% Creating equal loudness curve
fletcher_and_munson_40dB = fletcher_and_munson_40dB_curve_generator(frequency);

%% STEP 2
peaks = step2_interpolation(frequencyPeaksFiltered, ...
    magnitudeDecibels, frequency, fletcher_and_munson_40dB);
peaks = normalize_amplitudes(peaks, windowLength, beta);

%% STEP 3
Trajectories = step3_assign_peak_frequency_trajectories(peaks, MaximumPeakDeviation);

%% STEP 4
audioDataLength = length(audioData);
output = step4_resynthesize(Trajectories, fs, hopsize, audioDataLength);

audiowrite("output.wav",output,fs);

%% COMPARISONS
plot_results(audioData, output, fs, windowLength, beta, fftLength, overlapLength)

%% STEP 2
peaks = step2_interpolation_nofm40db(frequencyPeaksFiltered, ...
    magnitudeDecibels, frequency);
peaks = normalize_amplitudes(peaks, windowLength, beta);

%% STEP 3
Trajectories = step3_assign_peak_frequency_trajectories(peaks, MaximumPeakDeviation);

%% STEP 4
audioDataLength = length(audioData);
output = step4_resynthesize(Trajectories, fs, hopsize, audioDataLength);