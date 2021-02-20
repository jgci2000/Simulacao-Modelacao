clear all
close all
clc

A = 2;
fi = pi / 2;
%lambda = randi([1 6]);
lambda = 6;
x_1 = [0:0.01:12];
v = 0.1;
%t = 6;
x = 0:.001:16;


%Alínea a)

    x_2 = x_1 / lambda;
    y_1 = A * sin((2 * pi * x_2) + fi);
    plot(x_1, y_1, '-b', 'LineWidth', 1.5);
    %hold on
    %axis([x(i) x(i+12) -5 5]);

%Alínea d)
for t = 0:0.01:16
    
    %Alínea b)

    x_3 = (x_1 - v * t) / lambda;
    y_2 = A * sin((2 * pi * x_3) + fi);
    plot(x_1, y_2, '-y', 'LineWidth', 1.5);
    axis([0 12 -5 5]);
    hold on

    %Alínea c)

    x_4 = (x_1 + v * t) / lambda;
    y_3 = A * sin((2 * pi * x_4) + fi);
    plot(x_1, y_3, '-r', 'LineWidth', 1.5);

    y_4 = y_3 + y_2;
    plot(x_1, y_4, '-k', 'LineWidth', 3);
    pause(0.01);
    hold off
    
end








