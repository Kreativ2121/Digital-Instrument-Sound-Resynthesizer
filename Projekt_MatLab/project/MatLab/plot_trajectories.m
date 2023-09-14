clc;

for i = 1:size(Trajectories,2)
    x = [];
    y = [];
    for j = 1:4:size(Trajectories,1)
        if(isnan(Trajectories(j+1,i)) || isnan(Trajectories(j+3,i)) || (Trajectories(j+1,i)==0 && Trajectories(j+3,i)==0))
            continue;
        end
        x = [x, Trajectories(j+1,i)];
        y = [y, Trajectories(j+3,i)];
    end
    plot(x,y,'LineWidth',2);
    hold on;
end
title("Major Frequencies")
xlabel("Time (Trajectory no.)")
ylabel("Frequency (Hz)")
ylim([0 20000])
% xlim([0 105])
grid on