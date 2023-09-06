function plot_results(audioData, OutputAmp, fs, windowLength, beta, fftLength, overlapLength)
%compare_results Function that makes plots based on the input and output
%data of the algorythm.
%   This function makes three plots based on the data that was feeded to
%   algorythm and what came out of the algorythm. The first plot compares
%   2000 of points in a predefined range of the wavelet, to precisely check
%   the outcome data. The second plot shows both wavelets from beginning to
%   finish, to check for rare anomalies. The last plot compares STFT of the
%   original sound to the STFT of resynthesised sound.

    %COMPARE RESULTS
    f2 = figure('Name','Comparison','NumberTitle','off');
    f2.Position(1:2) = [650 850];
    subplot(2,1,1)
    plot(audioData(185000:187000))
    title("Original waveform")
    xlabel("Time (samples)")
    ylabel("Amplitude") 
    subplot(2,1,2)
    plot(OutputAmp(185000:187000))
    title("Resynthesized waveform")
    ylim([-1 1])
    xlabel("Time (samples)")
    ylabel("Amplitude")
    
    %COMPARE RESULTS OF WHOLE DATA
    f3 = figure('Name','Comparison of whole data','NumberTitle','off');
    f3.Position(1:2) = [50 55];
    f3.Position(3:4) = [2450 700];
    subplot(2,1,1)
    plot(audioData)
    title("Original waveform")
    xlabel("Time (samples)")
    ylabel("Amplitude") 
    subplot(2,1,2)
    plot(OutputAmp)
    title("Resynthesized waveform")
    ylim([-1 1])
    xlabel("Time (samples)")
    ylabel("Amplitude")
    
    %COMPARE ORIGINAL STFT WITH RESYNTHESIZED STFT
    f4 = figure('Name','STFT of resynthesized data','NumberTitle','off');
    f4.Position(1:2) = [1250 850];
    f4.Position(3:4) = [1250 420];
    subplot(2,1,1)
    stft(audioData, fs, Window=kaiser(windowLength,beta), ...
        FFTLength=fftLength, OverlapLength=overlapLength, ...
        FrequencyRange="onesided")
    title("STFT of original data")
    xlabel("Time (s)")
    ylabel("Frequency (kHz)")
    
    subplot(2,1,2)
    stft(OutputAmp, fs, Window=kaiser(windowLength,beta), ...
        FFTLength=fftLength, OverlapLength=overlapLength, ...
        FrequencyRange="onesided")
    title("STFT of resynthesized data")
    xlabel("Time (s)")
    ylabel("Frequency (kHz)")
end
