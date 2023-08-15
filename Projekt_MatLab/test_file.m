close all;
clear;
clc;

AmpSum = 0;
Output = [];

step = 1;
freq = 300;
% freq = [300,310,320,330,340,350,360,370,380,390,400];
for i = 1:30
    for j = 1:32
        if (i == 10)
            freq = freq+50;
        end
        AmpSum = 0;
        
        AmpSum = AmpSum + 0.8*cos((2*pi*freq*step)/44100);

        Output(step) = AmpSum;
        step = step + 1;

    end
end


plot(Output)
Output = Output';