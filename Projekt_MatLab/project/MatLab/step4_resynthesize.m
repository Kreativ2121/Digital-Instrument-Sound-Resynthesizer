function [output] = step4_resynthesize(trajectories, fs, hopsize, audioDataLength)
%step4_resynthesize Function that is responsible for resynthesizing sound.
%   This function takes data from spectral peak continuation and
%   reconstructs data by doing an additive synthesis per every time frame.
    output = [];
    stepcounter = 1;
    
    %Assign first synthesis frame
    nonZeroNonNaNtrajectories = numel(nonzeros(trajectories(1,:))) - sum(isnan(nonzeros(trajectories(1,:))));
    
    %Iterate over synth frames in first time frame
    for synframe = 1:hopsize
        ampSum = 0;
        % Iterate over peaks in a first synth frame
        for peak=1:nnz(trajectories(3,:))
            ampInst = trajectories(3,peak) + (trajectories(3,peak) ...
                -trajectories(3,peak))/nonZeroNonNaNtrajectories*synframe;
            freqInst = trajectories(4,peak);
            ampSum = ampSum + ampInst*sin((2*pi*freqInst*stepcounter)/fs);
        end
        output(stepcounter) = ampSum;
        stepcounter = stepcounter + 1;
    end
    
    %Iterate over time frames
    jumpMem = [];
    ampSumPrev = NaN;
    for tra = 2:size(trajectories,1)/4
        nonZeroTrajectories = numel(nonzeros(trajectories(tra*4-3,:)));
    
        %Iterate over synth frames
        synframesamount = floor(audioDataLength/(size(trajectories,1)/4));
        for synframe = 1:hopsize
    
            %Add peak sum calculated in the previous step
            ampSum = 0;
    
            %%Iterate over peaks
            for peak=1:nonZeroTrajectories

                % If trajectory is new (was non-existant in previous time frame)
                if((trajectories(((tra-2)*4)+3,peak)) == 0)
                    ampInst = 0 + (trajectories(((tra-1)*4)+3,peak))/synframesamount*(synframe-1);
                    freqInst = trajectories(((tra-1)*4)+4,peak);
                    ampSum = ampSum + ampInst*sin((2*pi*freqInst*stepcounter)/fs);
                    continue;
                % If trajectory died or is not born yet
                elseif(isnan(trajectories(((tra-1)*4)+3,peak)) || trajectories(((tra-1)*4)+3,peak)==0)
                    continue;
                % Traditional trajectory - not being born nor dying
                elseif(tra~=size(trajectories,1)/4)
                    ampInst = trajectories(((tra-2)*4)+3,peak) ...
                        +(trajectories(((tra-1)*4)+3,peak) ...
                        -trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                    freqInst = trajectories(((tra-2)*4)+4,peak)...
                        +(trajectories(((tra-1)*4)+4,peak)...
                        -trajectories(((tra-2)*4)+4,peak))/synframesamount*(synframe-1);
                    ampSum = ampSum + ampInst*sin((2*pi*freqInst*stepcounter)/fs);
                    continue;
                else
                    % Last time frame
                    ampInst = trajectories(((tra-2)*4)+3,peak) + (0 ...
                        -trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                    freqInst = trajectories(((tra-2)*4)+4,peak);
                    ampSum = ampSum + ampInst*sin((2*pi*freqInst*stepcounter)/fs);
                    continue;
                end
            end
    
            if(stepcounter >= 1000 && ampSumPrev<ampSum && abs(ampSum)<0.01)
                jumpMem = [jumpMem,stepcounter];
                stepcounter = 1;
                continue;
            end
    
            ampSumPrev = ampSum;
            output(sum(jumpMem)+stepcounter) = ampSum;
            stepcounter = stepcounter + 1;
        end
    end
    output = output';
    
    % FAILSAFE
    for i=1:size(output)
        if(output(i)>1)
           output(i) = 1.00;
        elseif(output(i)<-1)
            output(i) = -1.00;
        elseif(isnan(output(i)))
            output(i) = 0;
        end
    end
end

