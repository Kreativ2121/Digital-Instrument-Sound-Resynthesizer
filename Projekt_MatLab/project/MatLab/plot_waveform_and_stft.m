function plot_waveform_and_stft(audioData, fs, fftLength, windowLength, overlapLength, beta)
%plot_wavelet_and_stf Function that makes a plot of a few hundred first input audio
%values and makes a STFT plot with given variables.

    %% Plot wavelet
    xLimitation = [duration(0,0,0,0) duration(0,0,0,100)];
    
    f1 = figure('Name','STFT','NumberTitle','off');
    subplot(2,1,1)
    t = seconds(0:1/fs:(size(audioData,1)-1)/fs);
    f1.Position(1:2) = [50 850];
    plot(t,audioData)
    title("Waveform")
    xlabel("Time")
    ylabel("Amplitude")
    legend("Channel 1", "Channel 2")
    xlim(xLimitation)
    ylim([-1 1])
    
    %% Plot STFT
    subplot(2,1,2)
    stft(audioData, fs, Window=kaiser(windowLength,beta), ...
        FFTLength=fftLength, OverlapLength=overlapLength, ...
        FrequencyRange="onesided")
end

