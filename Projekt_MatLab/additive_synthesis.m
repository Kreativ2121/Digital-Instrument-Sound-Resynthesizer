clc;
clear;
xLimitation = [duration(0,0,0,0) duration(0,0,0,100)];
% xLimitation = "tight";

%% Ustawienia
MinPeak = 0.05;
MinPeakDist = 3;

%% Odczytaj plik wave,wygraj go na głośnikach i zapisz wave wyjściowy
% filetitle = "src/download/CantinaBand3.wav";
% filetitle = "src/records/kross/stereo/STRINGS_ACC_AMOL.wav";
% filetitle = "src/generated/mono/square440.wav";
% filetitle = "src/generated/mono/sine440.wav";
% filetitle = "src/generated/mono/tri440.wav";
filetitle = "src/generated/mono/saw440.wav";
[audioData,fs] = audioread(filetitle);
auInfo = audioinfo(filetitle);
% soundsc(audioData,fs)

% Zapis pliku wyjściowego (przepisanie pliku)
% audiowrite(filetitle,audioData,fs);

% Zapisanie do zmiennej tylko nazwy pliku ze ścieżki
file = strsplit(filetitle,'/');
file = file(end);

%% Wyplotuj przebieg sygnału z podziałem na kanały
t = seconds(0:1/fs:(size(audioData,1)-1)/fs);

subplot(3,1,1)
plot(t,audioData)
title("Przebieg")
xlabel("Czas")
ylabel("Amplituda")
legend("Kanał 1", "Kanał 2")
xlim(xLimitation)
%xlim("tight")
ylim([-1 1])

%% Transformacja Fouriera
subplot(3,1,2)

L = auInfo.TotalSamples;
Y = fft(audioData);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
%plot(f,P1)
% TODO Przed wyznaczaniem peaków można poddać funkcję filtracji podobnej do
% tych stosowanych w sygnałach EKG https://www.mathworks.com/help/signal/ref/findpeaks.html
findpeaks(P1,f,'MinPeakProminence',MinPeak,'MinPeakDistance',MinPeakDist)
[peaks,peaks2] = findpeaks(P1,f,'MinPeakProminence',MinPeak,'MinPeakDistance',MinPeakDist);
peaks(:,2) = peaks2';

title("FFT Spektrum")
xlabel("f (Hz)")
ylabel("|P1(f)|")

%% Synteza addytywna

subplot(3,1,3)
t = 0 : 1/fs : length(audioData)/fs;
t(end)=[];

signal = 0;
for i = 1:size(peaks,1)
    f1 = peaks(i,2); % per second
    T1 = 1/f1; % period, seconds
    amp1 = peaks(i,1); % amplitude
    
    signal1 = amp1 * sin(2*pi*t/T1);
    signal = signal + signal1;
    
end

plot(t, signal, 'b.-', 'LineWidth', 0.3, 'MarkerSize', 0.3);
xlim([0,0.1])
title("Zrekonstruowany sygnał")
xlabel("Czas")
ylabel("Amplituda")

%zapis do pliku
audiowrite(strcat("output/",file),signal,fs);

%% Odwrotna Transformacja Fouriera
% subplot(2,2,3)
% 
% Yi = ifft(Y, auInfo.TotalSamples);
% Timei = [];
% for c = 1:1:L
%     Timei(c) = c/fs;
% end
% 
% plot(Timei, Yi);
% title("IFFT Przebieg")
% xlabel("Czas (s)")
% ylabel("Amplituda")
% xlim([0,0.1])
% %xlim("tight")
% ylim([-1 1])

