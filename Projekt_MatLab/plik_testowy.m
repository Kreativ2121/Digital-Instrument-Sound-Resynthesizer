clc;
clear;

%% Odczytaj plik wave,wygraj go na głośnikach i zapisz wave wyjściowy
[audioData,fs] = audioread("src/download/CantinaBand3.wav");
soundsc(audioData,fs)
audiowrite("output/CantinaBand3.wav",audioData,fs);

%% Odczytywanie pliku frame-by-frame
% Tworzenie odczytywacza
fileReader = dsp.AudioFileReader("Filename","src/download/CantinaBand3.wav");
% Tworzenie writera, który będzie odtwarzał dźwięk
deviceWriter = audioDeviceWriter("SampleRate",fileReader.SampleRate);

while ~isDone(fileReader)
    audioData = fileReader();
    deviceWriter(audioData);
end

release(fileReader)
release(deviceWriter)