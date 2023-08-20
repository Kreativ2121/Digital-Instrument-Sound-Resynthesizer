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
    ampSumPrevPrev = NaN;
    
    for tra = 2:size(trajectories,1)/4
        nonZeroTrajectories = numel(nonzeros(trajectories(tra*4-3,:)));
    
        %Iterate over synth frames
        synframesamount = floor(audioDataLength/(size(trajectories,1)/4));
        for synframe = 1:hopsize
    
            %Add peak sum calculated in the previous step
            ampSum = 0;
    
            %%Iterate over peaks
            for peak=1:nonZeroTrajectories
    
                % Jeżeli ta trajektoria właśnie umarła
                if(isnan(trajectories(((tra-1)*4)+3,peak)) && ~isnan(trajectories(((tra-2)*4)+3,peak)))
                    ampInst = trajectories(((tra-2)*4)+3,peak) + (0 ...
                        -trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                    freqInst = trajectories(((tra-2)*4)+4,peak);
                    ampSum = ampSum + ampInst*sin((2*pi*freqInst*stepcounter)/fs);
                    continue;
    
                %Measure the instantaneous amplitude
                % Jeżeli nie istnieje już taka trajektoria, ale nie umarła w poprzednim framie.
                elseif(isnan(trajectories(((tra-1)*4)+3,peak)) || trajectories(((tra-1)*4)+3,peak)==0)
                    continue;
    
                % Jeżeli trajektoria jest nowa (nie istniała w poprzedniej próbce czasowej)
                elseif((trajectories(((tra-2)*4)+3,peak)) == 0)
                    ampInst = 0 + (trajectories(((tra-1)*4)+3,peak))/synframesamount*(synframe-1);
                    freqInst = trajectories(((tra-1)*4)+4,peak);
                    ampSum = ampSum + ampInst*sin((2*pi*freqInst*stepcounter)/fs);
    
                elseif(tra~=size(trajectories,1)/4)
                    if(isnan(trajectories(((tra)*4)+3,peak)))
                        % Tu wpisujemy dane trajektorii umierającej, ale jeszcze przed jej zgonem
                        ampInst = trajectories(((tra-2)*4)+3,peak) ...
                            +(trajectories(((tra-1)*4)+3,peak) ...
                            -trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                        freqInst = trajectories(((tra-2)*4)+4,peak) ...
                            +(trajectories(((tra-1)*4)+4,peak) ...
                            -trajectories(((tra-2)*4)+4,peak))/synframesamount*(synframe-1);
                        ampSum = ampSum + ampInst*sin((2*pi*freqInst*stepcounter)/fs);
                    else
                        % Zwykła trajektoria - nie rodząca się i nie umierająca
                        ampInst = trajectories(((tra-2)*4)+3,peak) ...
                            +(trajectories(((tra-1)*4)+3,peak) ...
                            -trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                        freqInst = trajectories(((tra-2)*4)+4,peak) + ...
                            (trajectories(((tra-1)*4)+4,peak)...
                            -trajectories(((tra-2)*4)+4,peak))/synframesamount*(synframe-1);
                        ampSum = ampSum + ampInst*sin((2*pi*freqInst*stepcounter)/fs);
                    end
                else
                    % Zwykła trajektoria - nie rodząca się i nie umierająca - na końcu wszystkich time frame'ów
                    ampInst = trajectories(((tra-2)*4)+3,peak) + ...
                        (trajectories(((tra-1)*4)+3,peak) - ...
                        trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                    freqInst = trajectories(((tra-2)*4)+4,peak);
                    ampSum = ampSum + ampInst*sin((2*pi*freqInst*stepcounter)/fs);
                end
            end
    
            % ampSumHist = [ampSumHist,ampSum];
            % if(stepcounter >= 10000 && ampSumPrev<ampSum && abs(ampSum)<0.005)
            % if(ampSumHist(end)<ampSum && abs(ampSum)<0.01)
            % if(ampSumPrev<ampSum && ampSumPrevPrev<ampSumPrev && abs(ampSum)<0.01)
            if(ampSumPrev<ampSum && ampSumPrevPrev<ampSumPrev && ...
                    abs(ampSum)<0.01 && stepcounter>= 100)
                jumpMem = [jumpMem,stepcounter];
                stepcounter = 1;
                continue;
            end
    
            ampSumPrevPrev = ampSumPrev;
            ampSumPrev = ampSum;
            output(sum(jumpMem)+stepcounter) = ampSum;
            stepcounter = stepcounter + 1;
        end
    end
    output = output';
    % output = output./3;
    
    % FAILSAFE
    for i=1:size(output)
        if(output(i)>1)
           output(i) = 1.00;
        elseif(output(i)<-1)
            output(i) = -1.00;
        end
    end
end

