clc;
clear;
close all;
% warning('on','verbose')
xLimitation = [duration(0,0,0,0) duration(0,0,0,100)];

%% Read Wave File
% filetitle = "src/generated/mono/square2000.wav";
% filetitle = "src/generated/mono/square440.wav";
% filetitle = "src/generated/mono/square689.wav";
% filetitle = "src/generated/mono/square2411.wav";
% filetitle = "src/generated/mono/sine440.wav";
% filetitle = "src/generated/mono/sine689.wav";
%filetitle = "src/generated/mono/saw689.wav";
% filetitle = "src/generated/mono/chirp440_2000.wav";
filetitle = "src/generated/mono/chirp2000_8000.wav";
% filetitle = "src/generated/mono/sine2000.wav";
% filetitle = "src/generated/mono/square2000_additivesynthesis.wav";
% filetitle = "src/download/CantinaBand3.wav";

[audioData,fs] = audioread(filetitle);
auInfo = audioinfo(filetitle);

% Zapisanie do zmiennej tylko nazwy pliku ze ścieżki
file = strsplit(filetitle,'/');
file = file(end);

%% Plot wavelet
t = seconds(0:1/fs:(size(audioData,1)-1)/fs);

f1 = figure('Name','STFT','NumberTitle','off');
f1.Position(1:2) = [50 850];
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

%% STEP 1 - STFT
% Default Hann128 window
% stft(audioData,fs)

% Kaiser window
fftlength = 128;
windowlength = fftlength;
%overlaplength = 96;
overlaplength = floor(0.75 * windowlength);
hopsize = fftlength - overlaplength;
beta = 6.0;
stft(audioData, fs, Window=kaiser(windowlength,beta), FFTLength=fftlength, OverlapLength=overlaplength, FrequencyRange="onesided") %Note
                                                                                                         %When this argument is set to "onesided", stft outputs the values in the positive 
                                                                                                         %Nyquist range and does not conserve the total power.
[magnitude,frequency,time] = stft(audioData,fs, Window=kaiser(windowlength,beta), FFTLength=fftlength, OverlapLength=overlaplength, FrequencyRange="onesided");

%% STEP 2 - CONVERSION TO POLAR COORDINATES
Ray = abs(magnitude);
PointAmplitude = deg2rad(angle(magnitude)); %TODO Czy ta linia ma sens?
% PointAmplitude = angle(magnitude);

%% STEP 3 - CHANGING MAGNITUDE SCALE TO DECIBELS
MagnitudeDecibels = 20*log10(Ray);

%% STEP 4 - FINDING PROMINENT SPECTRAL PEAKS
MinimumPeakHeightApprox = 20; %in dB
MinimumPeakHeightGlobal = -40; %in dB
MinimumPeakHeightLocal = -30; %in dB
FrequencyRangeLow = 20; %in hz
FrequencyRangeHigh = 16000; %in hz
AmplitudeRangeLow = -70; %in dB
AmplitudeRangeHigh = 0; %in dB

counter = 1;
FrequencyPeaks = [];

% Find peaks
% Loop on time frames
% TODO Sprawić by krawędziowe elementy mogły być wliczane jako piki
for j=1:size(MagnitudeDecibels,2)
    % Loop on frequency bins
    for i=2:size(MagnitudeDecibels,1)-1
        if(MagnitudeDecibels(i,j) >= MagnitudeDecibels(i-1,j) && MagnitudeDecibels(i,j) >= MagnitudeDecibels(i+1,j))
            FrequencyPeaks(1,counter) = MagnitudeDecibels(i,j);
            FrequencyPeaks(2,counter) = i; % No. of bin
            FrequencyPeaks(3,counter) = j; % No. of time frame
            FrequencyPeaks(7,counter) = PointAmplitude(i,j);
            counter = counter + 1;
        end
    end
end

% Find valleys for all peaks
PeakValleys = [];
for i=1:size(FrequencyPeaks,2)
    prevValue = FrequencyPeaks(1,i);
    for j = FrequencyPeaks(2,i):-1:1
        if(MagnitudeDecibels(j,FrequencyPeaks(3,i)) <= prevValue)
            prevValue = MagnitudeDecibels(j,FrequencyPeaks(3,i));
        else
            PeakValleys(i,1) = prevValue;
            break;
        end
    end
    PeakValleys(i,1) = prevValue;

    prevValue = FrequencyPeaks(1,i);
    for j = FrequencyPeaks(2,i):size(MagnitudeDecibels,1)
        if(MagnitudeDecibels(j,FrequencyPeaks(3,i)) <= prevValue)
            prevValue = MagnitudeDecibels(j,FrequencyPeaks(3,i));
        else
            PeakValleys(i,2) = prevValue;
            break;
        end
    end
    PeakValleys(i,2) = prevValue;
end
PeakValleys = PeakValleys';

% Assigning Peak Heights and unique index - usable for peak interpolation
for i=1:size(FrequencyPeaks,2)
    FrequencyPeaks(4,i) = FrequencyPeaks(1,i) - (PeakValleys(1,i) + PeakValleys(2,i))/2; % Peak heights in dB
    FrequencyPeaks(6,i) = i;
end

% % Adding magnitude data relative to maxdB (the "negative dB scale")
% Adding a row that will show values relative to max dB in whole sound
% maxdB = max(FrequencyPeaksRangeFiltered(4,:));
maxdB = max(FrequencyPeaks(1,:));
% maxdB = 0;

counter = 1;
for i=1:size(FrequencyPeaks,2)
    FrequencyPeaks(5,counter) = FrequencyPeaks(1,i)-maxdB;
    counter = counter+1;
end

% Adding the same row to Frequency Peaks (for Peak Interpolation) - needed?
% maxdB = max(FrequencyPeaks(4,:));
% counter = 1;
% for i=1:size(FrequencyPeaks,2)
%     FrequencyPeaks(5,counter) = FrequencyPeaks(1,i)-maxdB;
%     counter = counter+1;
% end

% Find maxAmplitudes in all Time Frames -> useful in discarding small peaks
 MaxAmplitudesInTimeFrames = zeros(1,size(MagnitudeDecibels,2));
for i=1:size(MagnitudeDecibels,2)
    MaxAmplitudesInTimeFrames(i) = max(MagnitudeDecibels(:,i)) - maxdB;
end

% Filtering Small Peaks
FrequencyPeaksHeightFiltered = [];
counter = 1;
for i=1:size(FrequencyPeaks,2)
    maxAmplitudeInTimeFrame = MaxAmplitudesInTimeFrames(FrequencyPeaks(3,i));
    MinimumPeakHeightLocalThreshold = maxAmplitudeInTimeFrame + MinimumPeakHeightLocal;

    % Filtrowanie zarówno na bazie obszaru w którym się znajduje amplituda,
    % jak i na bazie samej wielkości amplitudy pomniejszonej o doliny.
    if((FrequencyPeaks(5,i) >= MinimumPeakHeightGlobal || FrequencyPeaks(5,i) >= MinimumPeakHeightLocalThreshold) && (FrequencyPeaks(4,i) > MinimumPeakHeightApprox))
        FrequencyPeaksHeightFiltered(1,counter) = FrequencyPeaks(1,i);
        FrequencyPeaksHeightFiltered(2,counter) = FrequencyPeaks(2,i);
        FrequencyPeaksHeightFiltered(3,counter) = FrequencyPeaks(3,i);
        FrequencyPeaksHeightFiltered(4,counter) = FrequencyPeaks(4,i);
        FrequencyPeaksHeightFiltered(5,counter) = FrequencyPeaks(5,i);
        FrequencyPeaksHeightFiltered(6,counter) = FrequencyPeaks(6,i);
        FrequencyPeaksHeightFiltered(7,counter) = FrequencyPeaks(7,i);
        counter = counter+1;
    end
end

% Obtain frequency bins
% kl = FrequencyRangeLow*size(frequency,1)/fs;
% kh = FrequencyRangeHigh*size(frequency,1)/fs;

% Filtering peaks out of the audible frequency range
FrequencyPeaksRangeFiltered = [];
counter = 1;
for i=1:size(FrequencyPeaksHeightFiltered,2)
    if(frequency(FrequencyPeaksHeightFiltered(2,i)) >= FrequencyRangeLow && frequency(FrequencyPeaksHeightFiltered(2,i)) <= FrequencyRangeHigh)
        FrequencyPeaksRangeFiltered(1,counter) = FrequencyPeaksHeightFiltered(1,i);
        FrequencyPeaksRangeFiltered(2,counter) = FrequencyPeaksHeightFiltered(2,i);
        FrequencyPeaksRangeFiltered(3,counter) = FrequencyPeaksHeightFiltered(3,i);
        FrequencyPeaksRangeFiltered(4,counter) = FrequencyPeaksHeightFiltered(4,i);
        FrequencyPeaksRangeFiltered(5,counter) = FrequencyPeaksHeightFiltered(5,i);
        FrequencyPeaksRangeFiltered(6,counter) = FrequencyPeaksHeightFiltered(6,i);
        FrequencyPeaksRangeFiltered(7,counter) = FrequencyPeaksHeightFiltered(7,i);
        counter = counter+1;
    end
end

% Discard peaks with very low magnitude in general-dB-range
FrequencyPeaksdBFiltered = [];
counter = 1;
for i=1:size(FrequencyPeaksRangeFiltered,2)
    if(FrequencyPeaksRangeFiltered(4,i) >= AmplitudeRangeLow)
        FrequencyPeaksdBFiltered(1,counter) = FrequencyPeaksRangeFiltered(1,i);
        FrequencyPeaksdBFiltered(2,counter) = FrequencyPeaksRangeFiltered(2,i);
        FrequencyPeaksdBFiltered(3,counter) = FrequencyPeaksRangeFiltered(3,i);
        FrequencyPeaksdBFiltered(4,counter) = FrequencyPeaksRangeFiltered(4,i);
        FrequencyPeaksdBFiltered(5,counter) = FrequencyPeaksRangeFiltered(5,i);
        FrequencyPeaksdBFiltered(6,counter) = FrequencyPeaksRangeFiltered(6,i);
        FrequencyPeaksdBFiltered(7,counter) = FrequencyPeaksRangeFiltered(7,i);
        counter = counter+1;
    end
end

% Creating equal loudness curve
FletcherMundson40LikeLoudnessCurveX = double(.05 + double(4000./frequency));
FletcherMundson40LikeLoudnessCurveX = FletcherMundson40LikeLoudnessCurveX(ceil(length(FletcherMundson40LikeLoudnessCurveX)/2):end);
FletcherMundson40LikeLoudnessCurve = FletcherMundson40LikeLoudnessCurveX .* 10.^(-FletcherMundson40LikeLoudnessCurveX);

FletcherMundson40LikeLoudnessCurveFullX = double(.05 + double(4000./frequency));
FletcherMundson40LikeLoudnessCurveFull = FletcherMundson40LikeLoudnessCurveFullX .* 10.^(-FletcherMundson40LikeLoudnessCurveFullX);

for i=1:size(FrequencyPeaks,2)
    FrequencyPeaks(5,i) = FrequencyPeaks(5,i) - FletcherMundson40LikeLoudnessCurveFull(FrequencyPeaks(2,i));
    FrequencyPeaks(1,i) = FrequencyPeaks(1,i) - FletcherMundson40LikeLoudnessCurveFull(FrequencyPeaks(2,i));
end

%% STEP 5 - Peak detection and interpolation - for negative scaled dB magnitude
Peaks = [];
counter = 1;
for i=1:size(FrequencyPeaksdBFiltered,2)
    % Peaks(1,counter) = FrequencyPeaksPositive(1,i);
    % Peaks(2,counter) = FrequencyPeaksPositive(2,i);
    % Peaks(3,counter) = FrequencyPeaksPositive(3,i);
    % Peaks(4,counter) = FrequencyPeaksPositive(4,i);
    % Peaks(5,counter) = FrequencyPeaksPositive(5,i);
    % Peaks(6,counter) = FrequencyPeaksPositive(6,i);
    Peaks(1,counter) = FrequencyPeaksdBFiltered(1,i);
    Peaks(2,counter) = FrequencyPeaksdBFiltered(2,i);
    Peaks(3,counter) = FrequencyPeaksdBFiltered(3,i);
    Peaks(4,counter) = FrequencyPeaksdBFiltered(4,i);
    Peaks(5,counter) = FrequencyPeaksdBFiltered(5,i);
    Peaks(6,counter) = FrequencyPeaksdBFiltered(6,i);
    Peaks(13,counter) = FrequencyPeaksdBFiltered(7,i);

    % % TODO Sprawdzić czy dla alfy bety i gammy nie ujemnej wyniki będą
    % wychodzić poprawne.
    % ZMIANA 9.08 - usunięto "-maxdB" z równań na alfę, betę i gammę
    %%ALPHA
    if(FrequencyPeaksdBFiltered(2,i)-1 == 0)
        Peaks(7,counter) = (MagnitudeDecibels(FrequencyPeaksdBFiltered(2,i),FrequencyPeaksdBFiltered(3,i)) - FletcherMundson40LikeLoudnessCurveFull(FrequencyPeaksdBFiltered(2,i)));
        Peaks(14,counter) = PointAmplitude(FrequencyPeaksdBFiltered(2,i),FrequencyPeaksdBFiltered(3,i));
    else
        % Zakładamy, że valley jest takie samo dla elementów znajdujących się obok pików
        Peaks(7,counter) = (MagnitudeDecibels(FrequencyPeaksdBFiltered(2,i)-1,FrequencyPeaksdBFiltered(3,i)) - FletcherMundson40LikeLoudnessCurveFull(FrequencyPeaksdBFiltered(2,i)));
        Peaks(14,counter) = PointAmplitude(FrequencyPeaksdBFiltered(2,i)-1,FrequencyPeaksdBFiltered(3,i));
    end
    % angle(PointAmplitude(,))
    %%BETA
    Peaks(8,counter) = (MagnitudeDecibels(FrequencyPeaksdBFiltered(2,i),FrequencyPeaksdBFiltered(3,i)) - FletcherMundson40LikeLoudnessCurveFull(FrequencyPeaksdBFiltered(2,i)));
    Peaks(15,counter) = PointAmplitude(FrequencyPeaksdBFiltered(2,i),FrequencyPeaksdBFiltered(3,i));

    %%GAMMA
    if(FrequencyPeaksdBFiltered(2,i)+1 == length(FrequencyPeaksdBFiltered))
        Peaks(9,counter) = (MagnitudeDecibels(FrequencyPeaksdBFiltered(2,i),FrequencyPeaksdBFiltered(3,i)) - FletcherMundson40LikeLoudnessCurveFull(FrequencyPeaksdBFiltered(2,i)));
        Peaks(16,counter) = PointAmplitude(FrequencyPeaksdBFiltered(2,i),FrequencyPeaksdBFiltered(3,i));
    else
        Peaks(9,counter) = (MagnitudeDecibels(FrequencyPeaksdBFiltered(2,i)+1,FrequencyPeaksdBFiltered(3,i)) - FletcherMundson40LikeLoudnessCurveFull(FrequencyPeaksdBFiltered(2,i)));
        Peaks(16,counter) = PointAmplitude(FrequencyPeaksdBFiltered(2,i)+1,FrequencyPeaksdBFiltered(3,i));
    end

    %Assign parabola peak location - is it even necessary? - peak location
    %not given before?
    Peaks(10,counter) = ((Peaks(7,counter)-Peaks(9,counter))/(Peaks(7,counter)-2*Peaks(8,counter)+Peaks(9,counter)))/2;

    %Assign magnitude peak height
    Peaks(11,counter) = Peaks(8,counter) - ((Peaks(7,counter)-Peaks(9,counter))/4)*Peaks(10,counter);

    %Assign True Peak Location (in bins)
    Peaks(12,counter) = frequency(Peaks(2,counter)) + Peaks(10,counter);

    % % % PHASE CALCULATION - y(p)_theta
    Peaks(17,counter) = Peaks(15,counter) - ((Peaks(14,counter)-Peaks(16,counter))/4)*Peaks(10,counter);

    counter = counter+1;
end

% Normalization scale factor calculation for Kaiser window
N = size(MagnitudeDecibels,1);
I0_beta = besseli(0, beta); %TODO Sprawdzić poprawność wyniku.

Z = [];
for i = 1:128
    Z(i) = i-1;
end

% % SPRAWDZIĆ ZE WZOREM - CZYM JEST n?
W0_partial = real(besseli(0, beta .* sqrt(1-(2.*Z./(128-1))))./I0_beta);
W0 = sum(W0_partial);

% TODO W0 czy W0_B? - Które rozwiązanie lepsze?
W0_B = sum(kaiser(128,beta));
alpha = double(2/W0_B);

% Normalized amplitude
N_Amp = Peaks(11,:) * alpha; %Czy na pewno to ma być tak zmierzone? W końcu skala jest -70-0?

%% TODO SPRAWDZIĆ - Normalizacja zgodnie z punktem z ostatniego akapitu strony 47                          
Peaks(11,:) = N_Amp;

% % TEST - czemu wychodzi tak lepiej?
% Peaks(11,:) = Peaks(10,:);

%% STEP 6 - ASSIGNING PEAKS TO FREQUENCY TRAJECTORIES
tic

MaximumPeakDeviation = 500;
% MaximumPeakDeviation = 300; %Większa granica -> mniej trajektorii

PeaksMod = [];
PeaksMod(1,:) = Peaks(2,:);
PeaksMod(2,:) = Peaks(3,:);
PeaksMod(3,:) = Peaks(11,:);
PeaksMod(4,:) = Peaks(12,:);
% PeaksMod(5,:) = Peaks(13,:);
PeaksMod(5,:) = Peaks(17,:);
PeaksMod(6,:) = 0; % 0-unmatched 1-matched

% LowestNonNan = 1;
Trajectories = [];
FirstPeaksLoc = find(PeaksMod(2,:) == 1);
counter = 1;

% Assigning first time frame to trajectories
if(~isempty(FirstPeaksLoc))
    for i=1:length(FirstPeaksLoc)
        Trajectories(1,counter) = PeaksMod(1,FirstPeaksLoc(i));
        Trajectories(2,counter) = PeaksMod(2,FirstPeaksLoc(i));
        Trajectories(3,counter) = PeaksMod(3,FirstPeaksLoc(i));
        Trajectories(4,counter) = PeaksMod(4,FirstPeaksLoc(i));
        Trajectories(5,counter) = PeaksMod(5,FirstPeaksLoc(i));
        PeaksMod(6,i) = 1;
        counter = counter + 1;
    end
else
    Trajectories(1,counter) = NaN;
    Trajectories(2,counter) = NaN;
    Trajectories(3,counter) = NaN;
    Trajectories(4,counter) = NaN;
    Trajectories(5,counter) = NaN;
end

% Iterate over time frames
for i=2:max(PeaksMod(2,:))
    PeaksLoc = find(PeaksMod(2,:) == i);
    
    % Od razu poszerzamy macierz
    Trajectories(size(Trajectories,1)+5,1) = 0;
    % Trajectories(size(Trajectories,1)-4,1) = 0;
    % Trajectories(size(Trajectories,1)-3,1) = 0;
    % Trajectories(size(Trajectories,1)-2,1) = 0;
    % Trajectories(size(Trajectories,1)-1,1) = 0;

    NewPeaks = [];
    counter_in = 1;

    % FAILSAFE WHEN THERE ARE NO TRAJECTORIES IN A SINGLE FRAME
    if(isempty(PeaksLoc))
        NewPeaks(1, counter_in) = NaN;
        NewPeaks(2, counter_in) = NaN;
        NewPeaks(3, counter_in) = NaN;
        NewPeaks(4, counter_in) = NaN;
        NewPeaks(5, counter_in) = NaN;
        NewPeaks(6, counter_in) = NaN;
        Trajectories(4*(i-1)+1,1) = NaN;
        Trajectories(4*(i-1)+2,1) = NaN;
        Trajectories(4*(i-1)+3,1) = NaN;
        Trajectories(4*(i-1)+4,1) = NaN;
        Trajectories(4*(i-1)+5,1) = NaN;
        continue;
    end

    % Iterate over every peak in time frame
    for peak=PeaksLoc
        % Adding peaks from next time window to variable
        NewPeaks(1, counter_in) = peak;
        NewPeaks(2, counter_in) = PeaksMod(1,peak);
        NewPeaks(3, counter_in) = PeaksMod(2,peak);
        NewPeaks(4, counter_in) = PeaksMod(3,peak);
        NewPeaks(5, counter_in) = PeaksMod(4,peak);
        NewPeaks(6, counter_in) = PeaksMod(5,peak);
        PeaksMod(6, peak) = 1; % Mark as used.

        % Measuring distance from all previous peaks
        PeakLocCounter = 1;
        for j=7:7+length(Trajectories(size(Trajectories,1),:))-1
            NewPeaks(j, counter_in)=abs(PeaksMod(4,peak)-Trajectories(size(Trajectories,1)-6,PeakLocCounter)); %%Czy to -5 aby na pewno potrzebne?
            PeakLocCounter = PeakLocCounter + 1;
        end

        counter_in = counter_in + 1;
    end

    % If trajectory was previously killed, continue it with NaN
    for j = 1:size(Trajectories,2)
        if(isnan(Trajectories(size(Trajectories,1),j)))
            Trajectories(size(Trajectories,1)-4,j) = NaN;
            Trajectories(size(Trajectories,1)-3,j) = NaN;
            Trajectories(size(Trajectories,1)-2,j) = NaN;
            Trajectories(size(Trajectories,1)-1,j) = NaN;
            Trajectories(size(Trajectories,1),j) = NaN;
        end
    end
    
    % Choosing the smallest distance for every previous Trajectory
    PeakLocCounter = 1;
    TakenCounter = 1;
    
    condition = false;

    while condition ~= true
        closestVal = min(NewPeaks(7:end,1:end),[],"all");
        [minRow,minCol] = find(NewPeaks(7:end,1:end)==closestVal);
        
        if(closestVal > MaximumPeakDeviation)
            condition = true;

            continue;
        end
        
        Trajectories(5*(i-1)+1,minRow) = NewPeaks(2,minCol);
        Trajectories(5*(i-1)+2,minRow) = NewPeaks(3,minCol);
        Trajectories(5*(i-1)+3,minRow) = NewPeaks(4,minCol);
        Trajectories(5*(i-1)+4,minRow) = NewPeaks(5,minCol);
        Trajectories(5*(i-1)+5,minRow) = NewPeaks(6,minCol);
        % Discard the current peak
        for j=1:size(NewPeaks,1)
            NewPeaks(j,minCol) = NaN;
        end
        % Discard the trajectory associated with the current peak
        for j=1:size(NewPeaks,2)
            NewPeaks(6+minRow,j) = NaN;
        end
        PeakLocCounter = PeakLocCounter + 1;

        % If all trajectories have been set
        if(all(Trajectories(4*(i-1)+2,:)))
            condition = true;
            continue;
        end

        % If all peaks have been used
        if(all(all(isnan(NewPeaks(7:end,1:end)))))
            condition = true;
            continue;
        end
    end

    % Kill remaining trajectories
    for j = 1:size(Trajectories,2)
        if(Trajectories(size(Trajectories,1)-3,j) == 0)
            Trajectories(size(Trajectories,1)-4,j) = NaN;
            Trajectories(size(Trajectories,1)-3,j) = NaN;
            Trajectories(size(Trajectories,1)-2,j) = NaN;
            Trajectories(size(Trajectories,1)-1,j) = NaN;
            Trajectories(size(Trajectories,1),j) = NaN;
        end
    end

    % Create new trajectories from remaining peaks
    for j = 1:size(NewPeaks,2)
        if(~isnan(NewPeaks(1,j)))
            Trajectories(size(Trajectories,1)-4,size(Trajectories,2)+1) = NewPeaks(2,j);
            Trajectories(size(Trajectories,1)-3,size(Trajectories,2)) = NewPeaks(3,j);
            Trajectories(size(Trajectories,1)-2,size(Trajectories,2)) = NewPeaks(4,j);
            Trajectories(size(Trajectories,1)-1,size(Trajectories,2)) = NewPeaks(5,j);
            Trajectories(size(Trajectories,1),size(Trajectories,2)) = NewPeaks(6,j);
        end
    end
end
toc

%% STEP 7

OutputAmp = [];
stepcounter = 1;

%Assign first synthesis frame
NonZeroNonNaNTrajectories = numel(nonzeros(Trajectories(1,:))) - sum(isnan(nonzeros(Trajectories(1,:))));

%Iterate over synth frames in first time frame
for synframe = 1:hopsize
    AmpSum = 0;
    % Iterate over peaks in a first synth frame
    for peak=1:nnz(Trajectories(3,:))
        AmpInst = Trajectories(3,peak) + (Trajectories(3,peak) - Trajectories(3,peak))/NonZeroNonNaNTrajectories*synframe;
        % FreqInst = Trajectories(4,peak) + (Trajectories(4,peak) - Trajectories(4,peak))/NonZeroNonNaNTrajectories*synframe;
        FreqInst = Trajectories(4,peak);
        PhaseInst = Trajectories(5,peak);
        AmpSum = AmpSum + AmpInst*cos((2*pi*FreqInst*stepcounter)/fs);
    end
    OutputAmp(stepcounter) = AmpSum;
    stepcounter = stepcounter + 1;
end

%Iterate over time frames
AmpInst = 0;
AmpSum = 0;
AmpSumNext = [];
FreqInstNext = [];
NextTrajectory = [];
for tra = 2:size(Trajectories,1)/5
    NonZeroTrajectories = numel(nonzeros(Trajectories(tra*5-3,:)));
    NonZeroNonNaNTrajectories = numel(nonzeros(Trajectories(tra*5-3,:))) - sum(isnan(nonzeros(Trajectories(tra*5-3,:))));

    if(tra<(size(Trajectories,1)/5))
        NonZeroNonNaNTrajectoriesNext = numel(nonzeros(Trajectories((tra+1)*5-3,:))) - sum(isnan(nonzeros(Trajectories((tra+1)*5-3,:))));
    end

    %Iterate over synth frames
    synframesamount = floor(length(audioData)/(size(Trajectories,1)/5));
    for synframe = 1:hopsize

        %Add peak sum calculated in the previous step
        AmpSum = 0;

        %%Iterate over peaks
        for peak=1:NonZeroTrajectories
            
            % Jeżeli ta trajektoria właśnie umarła
            if(isnan(Trajectories(((tra-1)*5)+3,peak)) && ~isnan(Trajectories(((tra-2)*5)+3,peak)))
                AmpInst = Trajectories(((tra-2)*5)+3,peak) + (0 - Trajectories(((tra-2)*5)+3,peak))/synframesamount*synframe;

                % RÓŻNE PATENTY NA INTERPOLACJE:
                % FreqInst = Trajectories(((tra-2)*5)+4,peak)/synframesamount*synframe;
                % FreqInst = Trajectories(((tra-2)*5)+4,peak) + (0 - Trajectories(((tra-2)*5)+4,peak))/synframesamount*synframe;
                FreqInst = Trajectories(((tra-2)*5)+4,peak);
                
                AmpSum = AmpSum + AmpInst*cos((2*pi*FreqInst*stepcounter)/fs);
                continue;

            %Measure the instantaneous amplitude
            % Jeżeli nie istnieje już taka trajektoria, ale nie umarła w poprzednim framie.
            elseif(isnan(Trajectories(((tra-1)*5)+3,peak)) || Trajectories(((tra-1)*5)+3,peak)==0)
                continue;

            % Jeżeli trajektoria jest nowa (nie istniała w poprzedniej próbce czasowej)
            elseif((Trajectories(((tra-2)*5)+3,peak)) == 0)
                AmpInst = 0 + (Trajectories(((tra-1)*5)+3,peak))/synframesamount*synframe;
                % AmpInst = Trajectories(((tra-1)*5)+3,peak);
                % FreqInst = 0 + (Trajectories(((tra-1)*5)+4,peak))/synframesamount*synframe;
                FreqInst = Trajectories(((tra-1)*5)+4,peak); % Czy częstotliwość też mam interpolować, kiedy próbka wcześniej nie istniała?
                AmpSum = AmpSum + AmpInst*cos((2*pi*FreqInst*stepcounter)/fs);

            elseif(tra~=size(Trajectories,1)/5)
                if(isnan(Trajectories((tra*5)+3,peak)))
                    % Tu wpisujemy dane trajektorii umierającej, ale jeszcze przed jej zgonem
                    AmpInst = Trajectories(((tra-2)*5)+3,peak) + (Trajectories(((tra-1)*5)+3,peak) - Trajectories(((tra-2)*5)+3,peak))/synframesamount*synframe;
                    FreqInst = Trajectories(((tra-2)*5)+4,peak) + (Trajectories(((tra-1)*5)+4,peak) - Trajectories(((tra-2)*5)+4,peak))/synframesamount*synframe;
                    AmpSum = AmpSum + AmpInst*cos((2*pi*FreqInst*stepcounter)/fs);
                else
                    % Zwykła trajektoria - nie rodząca się i nie umierająca
                    AmpInst = Trajectories(((tra-2)*5)+3,peak) + (Trajectories(((tra-1)*5)+3,peak) - Trajectories(((tra-2)*5)+3,peak))/synframesamount*synframe;
                    % FreqInst = Trajectories(((tra-2)*5)+4,peak) + (Trajectories(((tra-1)*5)+4,peak) - Trajectories(((tra-2)*5)+4,peak))/synframesamount*synframe;
                    % FreqInst = Trajectories(((tra-1)*4)+4,peak);
                    % PhaseInst = Trajectories(5,peak);
                    % PhaseInst = phase_calculation_interpolation(Trajectories(((tra-2)*5)+5,peak),Trajectories(((tra-1)*5)+5,peak), Trajectories(((tra-2)*5)+4,peak), Trajectories(((tra-1)*5)+4,peak), synframe, synframesamount);
                    PhaseInst = phase_calculation_interpolation(Trajectories(((tra-2)*5)+5,peak),Trajectories(((tra-1)*5)+5,peak), 2*pi*Trajectories(((tra-2)*5)+4,peak), 2*pi*Trajectories(((tra-1)*5)+4,peak), synframe, synframesamount);
                    % PhaseInst = phase_calculation_interpolation(Trajectories(((tra-2)*5)+5,peak),Trajectories(((tra-1)*5)+5,peak), 2*pi*Trajectories(((tra-2)*5)+4,peak), 2*pi*Trajectories(((tra-1)*5)+4,peak), synframe, synframesamount);
                    AmpSum = AmpSum + AmpInst*cos(PhaseInst);
                end
            else
                % Zwykła trajektoria - nie rodząca się i nie umierająca - na końcu wszystkich time frame'ów
                AmpInst = Trajectories(((tra-2)*5)+3,peak) + (Trajectories(((tra-1)*5)+3,peak) - Trajectories(((tra-2)*5)+3,peak))/synframesamount*synframe;
                % FreqInst = Trajectories(((tra-2)*5)+4,peak) + (Trajectories(((tra-1)*5)+4,peak) - Trajectories(((tra-2)*5)+4,peak))/synframesamount*synframe;
                FreqInst = Trajectories(((tra-2)*5)+4,peak);
                AmpSum = AmpSum + AmpInst*cos((2*pi*FreqInst*stepcounter)/fs);
            end
        end
        OutputAmp(stepcounter) = AmpSum;
        stepcounter = stepcounter + 1;
    end
end
OutputAmp = OutputAmp';
% OutputAmp = OutputAmp./3;

% FAILSAFE
if(OutputAmp>1)
   OutputAmp = 1.00;
elseif(OutputAmp<-1)
    OutputAmp = -1.00;
end

%Zapisanie zresyntezowanego audio do pliku
% output = uint8(output);
audiowrite("output.wav",OutputAmp,fs);





% % % COMPARISONS

%COMPARE RESULTS
f2 = figure('Name','Comparison','NumberTitle','off');
f2.Position(1:2) = [650 850];
subplot(2,1,1)
plot(audioData(1:2000))
title("Przebieg oryginalny")
xlabel("Czas (próbki)")
ylabel("Amplituda") 
subplot(2,1,2)
plot(OutputAmp(1:2000))
title("Przebieg zresyntezowany")
ylim([-1 1])
xlabel("Czas (próbki)")
ylabel("Amplituda")

%COMPARE RESULTS OF WHOLE DATA
f3 = figure('Name','Comparison of whole data','NumberTitle','off');
f3.Position(1:2) = [50 55];
f3.Position(3:4) = [2450 700];
subplot(2,1,1)
plot(audioData)
title("Przebieg oryginalny")
xlabel("Czas (próbki)")
ylabel("Amplituda") 
subplot(2,1,2)
plot(OutputAmp)
title("Przebieg zresyntezowany")
ylim([-1 1])
xlabel("Czas (próbki)")
ylabel("Amplituda")

%COMPARE ORIGINAL STFT WITH RESYNTHESIZED STFT

f4 = figure('Name','STFT of resynthesized data','NumberTitle','off');
f4.Position(1:2) = [1250 850];
f4.Position(3:4) = [1250 420];
subplot(2,1,1)
stft(audioData, fs, Window=kaiser(windowlength,beta), FFTLength=fftlength, OverlapLength=overlaplength, FrequencyRange="onesided")
title("STFT of original data")
xlabel("Time (s)")
ylabel("Frequency (kHz)")

subplot(2,1,2)
stft(OutputAmp, fs, Window=kaiser(windowlength,beta), FFTLength=fftlength, OverlapLength=overlaplength, FrequencyRange="onesided")
title("STFT of resynthesized data")
xlabel("Time (s)")
ylabel("Frequency (kHz)")


%% Test 
% for i = 1:length(PeaksMod)
%     if(PeaksMod(2,i)~= i)
%         PeaksMod(2,i)
%     end
% end








