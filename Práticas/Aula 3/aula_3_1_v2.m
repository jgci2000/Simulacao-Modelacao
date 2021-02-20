close all
clear all

A = 2;
x = 0:0.001:12;
lambda = 8;
fi = pi / 6;
v = 0.1;

hold on

while lambda < 8.05
    while fi < (pi / 6) + 0.05
        for t = 0:0.1:6
            y_1 = A .* sin(2 .* pi .* (x ./ (lambda)) + fi);
            plot(x, y_1, '-b')
            xlim([0 12])
            x_2 = x - v * t;
            y_2 = A .* sin(2 .* pi .* (x_2 ./ lambda) + fi);
            plot(x, y_2, '-r')
            fi = fi + 0.1;
            lambda = lambda + 0.1;
            pause(0.1)
        end
    end
end
