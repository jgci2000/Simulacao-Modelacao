close all
clear all
clc

%Criação do triângulo
n = 3;
ang = 2 * pi / n;
res = zeros(2, n);
p0 = [3; 1];
R = [cos(ang), -sin(ang); sin(ang), cos(ang)];

for j = 1:n+1
    p0 = R * p0;
    res(:, j) = p0;
end

plot(res(1, :), res(2, :), '.y', 'MarkerSize', 35)
axis([-5 5 -5 5])
set(gca,'Color','k')
hold on

%Criação dos quadrados e rotação

for r = 0:0.5:6
    n = 4;
    ang = 2 * pi / n;
    res = zeros(2, n);
    p0 = [6; 6] - r;
    R = [cos(ang), -sin(ang); sin(ang), cos(ang)];

    for j = 1:n+1
        p0 = R * p0;
        res(:, j) = p0;
    end

    for i = 0:0.5:20
        plot(res(1, :) - i, res(2, :) - i, '*b', 'MarkerSize', 8)
        plot(res(1, :) + i, res(2, :) + i, '*b', 'MarkerSize', 8)
    end
end

%Criação do ponto (0,0)

p = plot(0, 0, '.g', 'MarkerSize', 20);
for i = 0:10
    set(p, 'Visible', 'off')
    pause(0.5)
    set(p, 'Visible', 'on')
    pause(0.5)
end
