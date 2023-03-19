clc;
clear;

%% Odczytaj plik wave,wygraj go na głośnikach i zapisz wave wyjściowy
filetitle = "src/download/CantinaBand3.wav";
[audioData,fs] = audioread(filetitle);
auInfo = audioinfo(filetitle);
% soundsc(audioData,fs)
audiowrite(filetitle,audioData,fs);

%% Wyplotuj przebieg sygnału z podziałem na kanały
t = seconds(0:1/fs:(size(audioData,1)-1)/fs);

subplot(3,1,1)
plot(t,audioData)
title("Przebieg")
xlabel("Czas")
ylabel("Amplituda")
legend("Kanał 1", "Kanał 2")
xlim("tight")
ylim([-1 1])

%% Wyplotuj obwiednię
[envMin,envMax,loc] = audioEnvelope(filetitle,NumPoints=2000);

subplot(3,1,2)
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

subplot(3,1,3)
loc = loc./fs;
plot(loc,envMax,loc,envMin)
title("Obwiednia 2",Interpreter="none")
xlabel("Czas")
ylabel("Amplituda")
xlim("tight")
ylim([-1 1])

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