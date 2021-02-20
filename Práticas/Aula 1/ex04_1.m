close all
clear all

A = 2;
x = 0:0.1:24;
lambda = 2;
while lambda < 6.1
    y = A .* cos(2 * pi * x / lambda);
    plot(x, y, 'lineWidth', 3);
    pause(0.3)
    lambda = lambda + 0.1;
end
