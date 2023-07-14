x = 0:0.01:1; % Zakres osi x
y = 4.2 * log10(0.1 * x + 1); % Funkcja aproksymująca

plot(x, y);
xlabel('Poziom dźwięku (x)');
ylabel('Głośność (dB)');
title('Aproksymacja krzywej Fletchera Mundsona');