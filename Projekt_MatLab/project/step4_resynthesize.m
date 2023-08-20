function [OutputAmp] = step4_resynthesize(Trajectories, fs, hopsize, audioDataLength)
%step4_resynthesize Summary of this function goes here
%   Detailed explanation goes here
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
            FreqInst = Trajectories(4,peak);
            AmpSum = AmpSum + AmpInst*sin((2*pi*FreqInst*stepcounter)/fs);
        end
        OutputAmp(stepcounter) = AmpSum;
        stepcounter = stepcounter + 1;
    end
    
    %Iterate over time frames
    jumpMem = [];
    
    AmpSumPrev = NaN;
    AmpSumPrevPrev = NaN;
    
    for tra = 2:size(Trajectories,1)/4
        NonZeroTrajectories = numel(nonzeros(Trajectories(tra*4-3,:)));
    
        %Iterate over synth frames
        synframesamount = floor(audioDataLength/(size(Trajectories,1)/4));
        for synframe = 1:hopsize
    
            %Add peak sum calculated in the previous step
            AmpSum = 0;
    
            %%Iterate over peaks
            for peak=1:NonZeroTrajectories
    
                % Jeżeli ta trajektoria właśnie umarła
                if(isnan(Trajectories(((tra-1)*4)+3,peak)) && ~isnan(Trajectories(((tra-2)*4)+3,peak)))
                    AmpInst = Trajectories(((tra-2)*4)+3,peak) + (0 - Trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                    FreqInst = Trajectories(((tra-2)*4)+4,peak);
                    AmpSum = AmpSum + AmpInst*sin((2*pi*FreqInst*stepcounter)/fs);
                    continue;
    
                %Measure the instantaneous amplitude
                % Jeżeli nie istnieje już taka trajektoria, ale nie umarła w poprzednim framie.
                elseif(isnan(Trajectories(((tra-1)*4)+3,peak)) || Trajectories(((tra-1)*4)+3,peak)==0)
                    continue;
    
                % Jeżeli trajektoria jest nowa (nie istniała w poprzedniej próbce czasowej)
                elseif((Trajectories(((tra-2)*4)+3,peak)) == 0)
                    AmpInst = 0 + (Trajectories(((tra-1)*4)+3,peak))/synframesamount*(synframe-1);
                    FreqInst = Trajectories(((tra-1)*4)+4,peak); % Czy częstotliwość też mam interpolować, kiedy próbka wcześniej nie istniała?
                    AmpSum = AmpSum + AmpInst*sin((2*pi*FreqInst*stepcounter)/fs);
    
                elseif(tra~=size(Trajectories,1)/4)
                    if(isnan(Trajectories(((tra)*4)+3,peak)))
                        % Tu wpisujemy dane trajektorii umierającej, ale jeszcze przed jej zgonem
                        AmpInst = Trajectories(((tra-2)*4)+3,peak) + (Trajectories(((tra-1)*4)+3,peak) - Trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                        FreqInst = Trajectories(((tra-2)*4)+4,peak) + (Trajectories(((tra-1)*4)+4,peak) - Trajectories(((tra-2)*4)+4,peak))/synframesamount*(synframe-1);
                        AmpSum = AmpSum + AmpInst*sin((2*pi*FreqInst*stepcounter)/fs);
                    else
                        % Zwykła trajektoria - nie rodząca się i nie umierająca
                        AmpInst = Trajectories(((tra-2)*4)+3,peak) + (Trajectories(((tra-1)*4)+3,peak) - Trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                        FreqInst = Trajectories(((tra-2)*4)+4,peak) + (Trajectories(((tra-1)*4)+4,peak) - Trajectories(((tra-2)*4)+4,peak))/synframesamount*(synframe-1);
                        AmpSum = AmpSum + AmpInst*sin((2*pi*FreqInst*stepcounter)/fs);
                    end
                else
                    % Zwykła trajektoria - nie rodząca się i nie umierająca - na końcu wszystkich time frame'ów
                    AmpInst = Trajectories(((tra-2)*4)+3,peak) + (Trajectories(((tra-1)*4)+3,peak) - Trajectories(((tra-2)*4)+3,peak))/synframesamount*(synframe-1);
                    FreqInst = Trajectories(((tra-2)*4)+4,peak);
                    AmpSum = AmpSum + AmpInst*sin((2*pi*FreqInst*stepcounter)/fs);
                end
            end
    
            % AmpSumHist = [AmpSumHist,AmpSum];
            % if(stepcounter >= 10000 && AmpSumPrev<AmpSum && abs(AmpSum)<0.005)
            % if(AmpSumHist(end)<AmpSum && abs(AmpSum)<0.01)
            % if(AmpSumPrev<AmpSum && AmpSumPrevPrev<AmpSumPrev && abs(AmpSum)<0.01)
            if(AmpSumPrev<AmpSum && AmpSumPrevPrev<AmpSumPrev && abs(AmpSum)<0.01 && stepcounter>= 100)
                jumpMem = [jumpMem,stepcounter];
                stepcounter = 1;
                continue;
            end
    
            AmpSumPrevPrev = AmpSumPrev;
            AmpSumPrev = AmpSum;
            OutputAmp(sum(jumpMem)+stepcounter) = AmpSum;
            stepcounter = stepcounter + 1;
        end
    end
    OutputAmp = OutputAmp';
    % OutputAmp = OutputAmp./3;
    
    % FAILSAFE
    for i=1:size(OutputAmp)
        if(OutputAmp(i)>1)
           OutputAmp(i) = 1.00;
        elseif(OutputAmp(i)<-1)
            OutputAmp(i) = -1.00;
        end
    end
end

