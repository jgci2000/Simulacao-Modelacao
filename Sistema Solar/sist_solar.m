clear all
close all
clc


for n = 1:1000
    
    ang = n * pi / 50;
    R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
    terra = [0; 0];
    terra = terra + 10;
    terra = R * terra;
    
    for i = 1:10
        hold off
        ang = i * pi / 25;
        R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
        lua = terra;
        lua = lua + 2;
        lua = R * lua;
        plot(lua(1, 1), lua(2, 1), 'ok', 'MarkerSize', 20, 'LineWidth', 2);
        sol = [0; 0];
        plot(sol(1, 1), sol(2, 1), 'oy', 'MarkerSize', 50, 'LineWidth', 5);
        hold on
        axis([-20 20 -20 20]);
        plot(terra(1, 1), terra(2, 1), 'ob', 'MarkerSize', 30, 'LineWidth', 3);
        pause(0.01);
    end
end