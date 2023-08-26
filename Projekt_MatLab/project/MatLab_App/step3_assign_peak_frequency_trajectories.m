function [trajectories] = step3_assign_peak_frequency_trajectories(peaks, maximumPeakDeviation)
%STEP3_ASSIGN_PEAK_FREQUENCY_trajectories A part of the algorithm that
%creates frequency trajectories based on the peaks feeded from the previous
%part of the algorythm.
%   This part of the code is responsible for the creation of "trajectories"
%   variable. This variable contains information about the frequency bin,
%   time frame, corrected frequency and amplitude of every major frequency
%   peak detected in the previous code. This variable is then responsible
%   for the resynthesis of the original sound.
    
    peaksMod = [];
    peaksMod(1,:) = peaks(2,:);
    peaksMod(2,:) = peaks(3,:);
    peaksMod(3,:) = peaks(11,:);
    peaksMod(4,:) = peaks(12,:);
    peaksMod(5,:) = 0; % 0-unmatched 1-matched
    
    % LowestNonNan = 1;
    trajectories = [];
    firstPeaksLoc = find(peaksMod(2,:) == 1);
    counter = 1;
    
    % Assigning first time frame to trajectories
    if(~isempty(firstPeaksLoc))
        for i=1:length(firstPeaksLoc)
            trajectories(1,counter) = peaksMod(1,firstPeaksLoc(i));
            trajectories(2,counter) = peaksMod(2,firstPeaksLoc(i));
            trajectories(3,counter) = peaksMod(3,firstPeaksLoc(i));
            trajectories(4,counter) = peaksMod(4,firstPeaksLoc(i));
            peaksMod(5,i) = 1;
            counter = counter + 1;
        end
    else
        trajectories(1,counter) = NaN;
        trajectories(2,counter) = NaN;
        trajectories(3,counter) = NaN;
        trajectories(4,counter) = NaN;
    end
    
    % Iterate over time frames
    for i=2:max(peaksMod(2,:))
        peaksLoc = find(peaksMod(2,:) == i);
        
        % Od razu poszerzamy macierz
        trajectories(size(trajectories,1)+4,1) = 0;
        trajectories(size(trajectories,1)-3,1) = 0;
        trajectories(size(trajectories,1)-2,1) = 0;
        trajectories(size(trajectories,1)-1,1) = 0;
    
        newPeaks = [];
        counter_in = 1;
    
        % FAILSAFE WHEN THERE ARE NO trajectories IN A SINGLE FRAME
        if(isempty(peaksLoc))
            newPeaks(1, counter_in) = NaN;
            newPeaks(2, counter_in) = NaN;
            newPeaks(3, counter_in) = NaN;
            newPeaks(4, counter_in) = NaN;
            newPeaks(5, counter_in) = NaN;
            trajectories(4*(i-1)+1,1) = NaN;
            trajectories(4*(i-1)+2,1) = NaN;
            trajectories(4*(i-1)+3,1) = NaN;
            trajectories(4*(i-1)+4,1) = NaN;
            continue;
        end
    
        % Iterate over every peak in time frame
        for peak=peaksLoc
            % Adding peaks from next time window to variable
            newPeaks(1, counter_in) = peak;
            newPeaks(2, counter_in) = peaksMod(1,peak);
            newPeaks(3, counter_in) = peaksMod(2,peak);
            newPeaks(4, counter_in) = peaksMod(3,peak);
            newPeaks(5, counter_in) = peaksMod(4,peak);
            peaksMod(5, peak) = 1; % Mark as used.
    
            peakLocCounter = 1;
            for j=6:6+length(trajectories(size(trajectories,1),:))-1
                newPeaks(j, counter_in)=abs(peaksMod(4,peak)- ...
                    trajectories(size(trajectories,1)-4,peakLocCounter));
                peakLocCounter = peakLocCounter + 1;
            end
    
            counter_in = counter_in + 1;
        end
    
        % If trajectory was previously killed, continue it with NaN
        for j = 1:size(trajectories,2)
            if(isnan(trajectories(size(trajectories,1),j)))
                trajectories(size(trajectories,1)-3,j) = NaN;
                trajectories(size(trajectories,1)-3,j) = NaN;
                trajectories(size(trajectories,1)-3,j) = NaN;
                trajectories(size(trajectories,1)-3,j) = NaN;
            end
        end
        
        % Choosing the smallest distance for every previous Trajectory
        peakLocCounter = 1;
        
        condition = false;
    
        while condition ~= true
            closestVal = min(newPeaks(6:end,1:end),[],"all");
            [minRow,minCol] = find(newPeaks(6:end,1:end)==closestVal);
            
            if(closestVal > maximumPeakDeviation)
                condition = true;
    
                continue;
            end
            
            trajectories(4*(i-1)+1,minRow) = newPeaks(2,minCol);
            trajectories(4*(i-1)+2,minRow) = newPeaks(3,minCol);
            trajectories(4*(i-1)+3,minRow) = newPeaks(4,minCol);
            trajectories(4*(i-1)+4,minRow) = newPeaks(5,minCol);
            
            % Discard the current peak
            for j=1:size(newPeaks,1)
                newPeaks(j,minCol) = NaN;
            end
            
            % Discard the trajectory associated with the current peak
            for j=1:size(newPeaks,2)
                newPeaks(5+minRow,j) = NaN;
            end
            peakLocCounter = peakLocCounter + 1;
    
            % If all trajectories have been set
            if(all(trajectories(4*(i-1)+2,:)))
                condition = true;
                continue;
            end
    
            % If all peaks have been used
            if(all(all(isnan(newPeaks(6:end,1:end)))))
                condition = true;
                continue;
            end
        end
    
        % Kill remaining trajectories
        for j = 1:size(trajectories,2)
            if(trajectories(size(trajectories,1)-3,j) == 0)
                trajectories(size(trajectories,1)-3,j) = NaN;
                trajectories(size(trajectories,1)-2,j) = NaN;
                trajectories(size(trajectories,1)-1,j) = NaN;
                trajectories(size(trajectories,1),j) = NaN;
            end
        end
    
        % Create new trajectories from remaining peaks
        for j = 1:size(newPeaks,2)
            if(~isnan(newPeaks(1,j)))
                trajectories(size(trajectories,1)-3,size(trajectories,2)+1) = newPeaks(2,j);
                trajectories(size(trajectories,1)-2,size(trajectories,2)) = newPeaks(3,j);
                trajectories(size(trajectories,1)-1,size(trajectories,2)) = newPeaks(4,j);
                trajectories(size(trajectories,1),size(trajectories,2)) = newPeaks(5,j);
            end
        end
    end
end

