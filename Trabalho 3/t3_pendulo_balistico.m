clear all
close all
clc


%% Constants

% Bullet
m1 = 0.100;
v1 = 250;

% Pendulum
m2 = 15;
v2 = 0;
L = 4;
g = 9.81;
theta0 = acos(-((v1 / (1 + (m2 / m1))) ^2 / (2 * g * L)) + 1);

%Constants

%% ODE Solver - Euler-Cromer Method

theta(1) = theta0;
dTheta(1) = 0;
tf = 2;
dt = 0.01;
t = 0:dt:tf;

for i = 1:length(t) - 1
    dTheta(i + 1) = dTheta(i) - (g / L )* theta(i) * dt;
    theta(i + 1) = theta(i) + dTheta(i + 1) * dt;
end

% plot(t, theta);

%% Equações de movimento

for i = 1:length(t)
    x = L * sin(theta(i));
    y = L * cos(theta(i));
    
    plot(x, y, '.k', 'MarkerSize', 100)
    hold on
    line([x 0], [y 0], 'LineWidth', 2)
    axis([-10 10 -10 10])
    axis equal
    hold off
end








