clc;
clear;

%% Odczytaj plik wave,wygraj go na głośnikach i zapisz wave wyjściowy
% filetitle = "src/download/CantinaBand3.wav";
% filetitle = "src/records/kross/stereo/STRINGS_ACC_AMOL.wav";
filetitle = "src/generated/mono/STRINGS_ACC_AMOL.wav";
[audioData,fs] = audioread(filetitle);
auInfo = audioinfo(filetitle);
soundsc(audioData,fs)
audiowrite(filetitle,audioData,fs);

%% Odczytywanie pliku frame-by-frame
% % Tworzenie odczytywacza
% fileReader = dsp.AudioFileReader("Filename","src/download/CantinaBand3.wav");
% % Tworzenie writera, który będzie odtwarzał dźwięk
% deviceWriter = audioDeviceWriter("SampleRate",fileReader.SampleRate);
% 
% while ~isDone(fileReader)
%     audioData = fileReader();
%     % deviceWriter(audioData);
% end
% 
% release(fileReader)
% release(deviceWriter)

%% Wyplotuj przebieg sygnału z podziałem na kanały
t = seconds(0:1/fs:(size(audioData,1)-1)/fs);

subplot(3,3,4)
plot(t,audioData)
title("Przebieg")
xlabel("Czas")
ylabel("Amplituda")
legend("Kanał 1", "Kanał 2")
xlim("tight")
ylim([-1 1])

%% Wyplotuj obwiednię
[envMin,envMax,loc] = audioEnvelope(filetitle,NumPoints=2000);

subplot(3,3,1)
nChans = size(envMin,2);
envbars = [shiftdim(envMin,-1);
    shiftdim(envMax,-1);
    shiftdim(NaN(size(envMin)),-1)];
ybars = reshape(envbars,[],nChans);
t = seconds(loc/auInfo.SampleRate);
tbars = reshape(repmat(t,3,1),[],1);
plot(tbars,ybars);
title("Obwiednia",Interpreter="none")
xlabel("Czas")
ylabel("Amplituda")
xlim("tight")
ylim([-1 1])

subplot(3,3,7)
loc = loc./fs;
plot(loc,envMax,loc,envMin)
title("Obwiednia 2",Interpreter="none")
xlabel("Czas")
ylabel("Amplituda")
xlim("tight")
ylim([-1 1])

%% Transformacja Fouriera
subplot(3,3,2)

L = auInfo.TotalSamples;
Y = fft(audioData);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
plot(f,P1)
title("FFT Spektrum")
xlabel("f (Hz)")
ylabel("|P1(f)|")

%% Odwrotna Transformacja Fouriera
subplot(3,3,5)

Yi = ifft(Y, auInfo.TotalSamples);
Timei = [];
for c = 1:1:L
    Timei(c) = c/fs;
end

plot(Timei, Yi);
title("IFFT Przebieg")
xlabel("Czas (s)")
ylabel("Amplituda")
xlim("tight")
ylim([-1 1])

%% STFT
subplot(3,3,3)

L = auInfo.TotalSamples;
Y = stft(audioData);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
plot(f,P1)
title("STFT Spektrum")
xlabel("f (Hz)")
ylabel("|P1(f)|")

%% ISTFT
subplot(3,3,6)

Yi = istft(Y, auInfo.TotalSamples);
Timei = [];
for c = 1:1:L
    Timei(c) = c/fs;
end

% Padding - co by istft zadziałało
while(length(Timei)>length(Yi))
    Yi(end+1) = 0;
end


plot(Timei, Yi);
title("ISTFT Przebieg")
xlabel("Czas (s)")
ylabel("Amplituda")
xlim("tight")
ylim([-1 1])
















