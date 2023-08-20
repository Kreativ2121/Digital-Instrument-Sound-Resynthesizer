function [Peaks] = step2_peak_detection_and_interpolation(FrequencyPeaksdBFiltered, magnitude, magnitudeDecibels, frequency, fletcher_and_mundson_40dB)
%step2_peak_detection_and_interpolation Summary of this function goes here
%   Detailed explanation goes here
    Peaks = [];
    distanceBetweenFreqBins = frequency(2)-frequency(1);
    counter = 1;
    for i=1:size(FrequencyPeaksdBFiltered,2)
        Peaks(1,counter) = FrequencyPeaksdBFiltered(1,i);
        Peaks(2,counter) = FrequencyPeaksdBFiltered(2,i);
        Peaks(3,counter) = FrequencyPeaksdBFiltered(3,i);
        Peaks(4,counter) = FrequencyPeaksdBFiltered(4,i);
        Peaks(5,counter) = FrequencyPeaksdBFiltered(5,i);
        Peaks(6,counter) = FrequencyPeaksdBFiltered(6,i);
    
        %%ALPHA
        if(FrequencyPeaksdBFiltered(2,i)-1 == 0)
            Peaks(7,counter) = (magnitudeDecibels(FrequencyPeaksdBFiltered(2,i),FrequencyPeaksdBFiltered(3,i)) - fletcher_and_mundson_40dB(FrequencyPeaksdBFiltered(2,i)));
        else
            % Zakładamy, że valley jest takie samo dla elementów znajdujących się obok pików
            Peaks(7,counter) = (magnitudeDecibels(FrequencyPeaksdBFiltered(2,i)-1,FrequencyPeaksdBFiltered(3,i)) - fletcher_and_mundson_40dB(FrequencyPeaksdBFiltered(2,i)));
        end
    
        %%BETA
        Peaks(8,counter) = (magnitudeDecibels(FrequencyPeaksdBFiltered(2,i),FrequencyPeaksdBFiltered(3,i)) - fletcher_and_mundson_40dB(FrequencyPeaksdBFiltered(2,i)));
    
        %%GAMMA
        if(FrequencyPeaksdBFiltered(2,i)+1 == length(FrequencyPeaksdBFiltered))
            Peaks(9,counter) = (magnitudeDecibels(FrequencyPeaksdBFiltered(2,i),FrequencyPeaksdBFiltered(3,i)) - fletcher_and_mundson_40dB(FrequencyPeaksdBFiltered(2,i)));
        else
            Peaks(9,counter) = (magnitudeDecibels(FrequencyPeaksdBFiltered(2,i)+1,FrequencyPeaksdBFiltered(3,i)) - fletcher_and_mundson_40dB(FrequencyPeaksdBFiltered(2,i)));
        end
    
        %Assign parabola peak location
        Peaks(10,counter) = ((Peaks(7,counter)-Peaks(9,counter))/(Peaks(7,counter)-2*Peaks(8,counter)+Peaks(9,counter)))/2;
    
        %Assign magnitude peak height
        Peaks(11,counter) = Peaks(8,counter) - ((Peaks(7,counter)-Peaks(9,counter))/4)*Peaks(10,counter);
    
        %Assign True Peak Location (in bins)
        Peaks(12,counter) = frequency(Peaks(2,counter)) + Peaks(10,counter)*distanceBetweenFreqBins;
    
        cplx = magnitude(Peaks(2,counter),Peaks(3,counter));
    
        if(Peaks(3,counter)==1)
            cplx_last = magnitude(Peaks(2,counter),Peaks(3,counter));
        else
            cplx_last = magnitude(Peaks(2,counter),Peaks(3,counter)-1);
        end
    
        if(Peaks(3,counter)==size(magnitude,2))
            cplx_next = magnitude(Peaks(2,counter),Peaks(3,counter));
        else
            cplx_next = magnitude(Peaks(2,counter),Peaks(3,counter)+1);
        end
    
        real_yp = real(20*log10(cplx)) - (1/4)*(real(20*log10(cplx_last))-real(20*log10(cplx_next)))*Peaks(10,counter);
        imag_yp = imag(20*log10(cplx)) - (1/4)*(imag(20*log10(cplx_last))-imag(20*log10(cplx_next)))*Peaks(10,counter);
        
        new_cplx = complex(real_yp,imag_yp);
    
        counter = counter+1;
    end
end

