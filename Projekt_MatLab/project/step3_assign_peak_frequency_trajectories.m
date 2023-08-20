function [Trajectories] = step3_assign_peak_frequency_trajectories(Peaks, MaximumPeakDeviation)
%STEP3_ASSIGN_PEAK_FREQUENCY_TRAJECTORIES Summary of this function goes here
%   Detailed explanation goes here
    
    PeaksMod = [];
    PeaksMod(1,:) = Peaks(2,:);
    PeaksMod(2,:) = Peaks(3,:);
    PeaksMod(3,:) = Peaks(11,:);
    PeaksMod(4,:) = Peaks(12,:);
    PeaksMod(5,:) = 0; % 0-unmatched 1-matched
    
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
            PeaksMod(5,i) = 1;
            counter = counter + 1;
        end
    else
        Trajectories(1,counter) = NaN;
        Trajectories(2,counter) = NaN;
        Trajectories(3,counter) = NaN;
        Trajectories(4,counter) = NaN;
    end
    
    % Iterate over time frames
    for i=2:max(PeaksMod(2,:))
        PeaksLoc = find(PeaksMod(2,:) == i);
        
        % Od razu poszerzamy macierz
        Trajectories(size(Trajectories,1)+4,1) = 0;
        Trajectories(size(Trajectories,1)-3,1) = 0;
        Trajectories(size(Trajectories,1)-2,1) = 0;
        Trajectories(size(Trajectories,1)-1,1) = 0;
    
        NewPeaks = [];
        counter_in = 1;
    
        % FAILSAFE WHEN THERE ARE NO TRAJECTORIES IN A SINGLE FRAME
        if(isempty(PeaksLoc))
            NewPeaks(1, counter_in) = NaN;
            NewPeaks(2, counter_in) = NaN;
            NewPeaks(3, counter_in) = NaN;
            NewPeaks(4, counter_in) = NaN;
            NewPeaks(5, counter_in) = NaN;
            Trajectories(4*(i-1)+1,1) = NaN;
            Trajectories(4*(i-1)+2,1) = NaN;
            Trajectories(4*(i-1)+3,1) = NaN;
            Trajectories(4*(i-1)+4,1) = NaN;
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
            PeaksMod(5, peak) = 1; % Mark as used.
    
            PeakLocCounter = 1;
            for j=6:6+length(Trajectories(size(Trajectories,1),:))-1
                % NewPeaks(j, counter_in)=abs(PeaksMod(4,peak)-Trajectories(4+(4*(i-2)),PeakLocCounter));
                NewPeaks(j, counter_in)=abs(PeaksMod(4,peak)-Trajectories(size(Trajectories,1)-4,PeakLocCounter));
                PeakLocCounter = PeakLocCounter + 1;
            end
    
            counter_in = counter_in + 1;
        end
    
        % If trajectory was previously killed, continue it with NaN
        for j = 1:size(Trajectories,2)
            if(isnan(Trajectories(size(Trajectories,1),j)))
                Trajectories(size(Trajectories,1)-3,j) = NaN;
                Trajectories(size(Trajectories,1)-3,j) = NaN;
                Trajectories(size(Trajectories,1)-3,j) = NaN;
                Trajectories(size(Trajectories,1)-3,j) = NaN;
            end
        end
        
        % Choosing the smallest distance for every previous Trajectory
        PeakLocCounter = 1;
        TakenCounter = 1;
        
        condition = false;
    
        while condition ~= true
            closestVal = min(NewPeaks(6:end,1:end),[],"all");
            [minRow,minCol] = find(NewPeaks(6:end,1:end)==closestVal);
            
            if(closestVal > MaximumPeakDeviation)
                condition = true;
    
                continue;
            end
            
            Trajectories(4*(i-1)+1,minRow) = NewPeaks(2,minCol);
            Trajectories(4*(i-1)+2,minRow) = NewPeaks(3,minCol);
            Trajectories(4*(i-1)+3,minRow) = NewPeaks(4,minCol);
            Trajectories(4*(i-1)+4,minRow) = NewPeaks(5,minCol);
            
            % Discard the current peak
            for j=1:size(NewPeaks,1)
                NewPeaks(j,minCol) = NaN;
            end
            
            % Discard the trajectory associated with the current peak
            for j=1:size(NewPeaks,2)
                NewPeaks(5+minRow,j) = NaN;
            end
            PeakLocCounter = PeakLocCounter + 1;
    
            % If all trajectories have been set
            if(all(Trajectories(4*(i-1)+2,:)))
                condition = true;
                continue;
            end
    
            % If all peaks have been used
            if(all(all(isnan(NewPeaks(6:end,1:end)))))
                condition = true;
                continue;
            end
        end
    
        % Kill remaining trajectories
        for j = 1:size(Trajectories,2)
            if(Trajectories(size(Trajectories,1)-3,j) == 0)
                Trajectories(size(Trajectories,1)-3,j) = NaN;
                Trajectories(size(Trajectories,1)-2,j) = NaN;
                Trajectories(size(Trajectories,1)-1,j) = NaN;
                Trajectories(size(Trajectories,1),j) = NaN;
            end
        end
    
        % Create new trajectories from remaining peaks
        for j = 1:size(NewPeaks,2)
            if(~isnan(NewPeaks(1,j)))
                Trajectories(size(Trajectories,1)-3,size(Trajectories,2)+1) = NewPeaks(2,j);
                Trajectories(size(Trajectories,1)-2,size(Trajectories,2)) = NewPeaks(3,j);
                Trajectories(size(Trajectories,1)-1,size(Trajectories,2)) = NewPeaks(4,j);
                Trajectories(size(Trajectories,1),size(Trajectories,2)) = NewPeaks(5,j);
            end
        end
    end
end

