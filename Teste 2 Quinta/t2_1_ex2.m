close all
clear all
clc


T1 = [0, 230; 1, 300; 2.5, 335; 3, 340];
T2 = [0.5, 300; 1.5, 350; 2, 390; 3, 400];

t = 0:0.5:3;

inter1 = spline(T1(:, 1), T1(:, 2), t);
inter2 = spline(T2(:, 1), T2(:, 2), t);

q = inter1 ./ inter2;

plot(t, q, '-b')

% plot(t, interp1, '-b')
% plot(t, interp2, '-r')
