clear all
close all
clc

% x0 = 0.5;
% R = 1;
% x = R * x0 * (1 - x0);
% mat = zeros(1, 102);
% mat(1, 1) = x0;
% mat(1, 2) = x;

for j = 1:0.01:4
    R = j;
    x0 = 0.5;
    x = R * x0 * (1 - x0); 
    mat = zeros(1, 1002);
    mat(1, 1) = x0;
    mat(1, 2) = x;
    for i = 1:1000
        x_ant = x;
        x = x_ant * R * (1 - x_ant);
        mat(1, i + 2) = x;
    end
    figure(1)
    plot(mat(1, :), '-r');
    axis([0 1002 0 1])
    %pause(0.05);
    figure(2)
    plot(R*ones(51), mat(1, end-50:end), '.b', 'MarkerSize', 1);
    hold on
end
%plot(mat(1, :), '-r');
%axis([0 (i + 2) 0 100])


%plot(R, mat(1, [end-20 end]), '*b');