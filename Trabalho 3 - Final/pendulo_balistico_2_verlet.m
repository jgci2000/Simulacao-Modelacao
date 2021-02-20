%% Constants

clear all
close all
clc

% Pendulum
m2 = 15;
v2 = 0;
L = 10;
x = 0; y = -L;

% Bullet
m1 = 0.3;
v1 = 250;
xb = -10;
yb = -L;

% Constants
g = 9.81;
theta0 = acos(-((v1 / (1 + (m2 / m1))) ^2 / (2 * g * L)) + 1);

b = 0.15;
M = m1 + m2;

% Plot Constants

tempo = 40;
nPoints = round(tempo * 1001 / 40);
dt = 0.04;
omega = zeros(nPoints, 1);
theta = zeros(nPoints, 1);
t = zeros(nPoints, 1);
theta(1) = theta0;
theta(2) = theta0 - 0.0001;

%% ODE Solver - Verlet Method
if b == 0
    for i = 1:nPoints - 2
        theta(i + 2) = 2 * theta(i + 1) - theta(i) - (g / L) * theta(i + 1) * dt ^2;
        omega(i + 1) = (theta(i + 2) - theta(i)) / (2 * dt);
        
        t(i + 2) = t(i + 1) + dt;
    end
elseif not(b == 0)
    for i = 1:nPoints - 2
        theta(i + 2) = 2 * theta(i + 1) - theta(i) - dt ^2 * ((g / L) * theta(i + 1) + b * omega(i));  % Last Expression is for the Damping Strength
        omega(i + 1) = (theta(i + 2) - theta(i)) / (2 * dt);
        
        t(i + 2) = t(i + 1) + dt;
    end
end

% 1001 points equals to 40 secunds
% plot(t, theta, 'r')

%% Animation Plots

figure(5)
% Bullet Impact
while xb < 0
    tf = (x - xb) / v1;
    dt2 = 0.001;
    
    for t1 = 0:dt2:tf
        xb = xb + v1 * dt2;
        
        plot(xb, yb, '.r', 'MarkerSize', 20)
        title("Simulação do Pendulo")
        hold on
        plot(x, y, '.r', 'MarkerSize', 100)
        line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
        axis([-10 10 -L-2 0])
        drawnow update
        hold off
    end
end

% From (0, 0) to initial position
for dTheta = 0:0.02:theta0
    x = L * sin(dTheta);
    y = - L * cos(dTheta);
    
    plot(x, y, '.r', 'MarkerSize', 100)
    title("Simulação do Pendulo")
    hold on
    line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
    axis([-10 10 -L-2 0])
    drawnow update
    hold off
end

% Simple Pendulum Animation
for i = 1:nPoints
    x = L * sin(theta(i));
    y = - L * cos(theta(i));
    
    plot(x, y, '.r', 'MarkerSize', 100)
    title("Simulação do Pendulo")
    line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
    axis([-10 10 -L-2 0])
    drawnow update
    hold off
end

%% Theta functino of time
figure(1)
% subplot(4, 1, 1)
plot(t, theta)
title("Variação de theta ao longo do tempo")
xlabel("t / s")
ylabel("theta / rad")
%%
% %% Energys

% %% Kinetic
% figure(2)
% subplot(4, 1, 2)
K = 0.5 * M * (omega .* L) .^2;
plot(t, K)
title("Energia Cinética ao longo do tempo")
xlabel("t / s")
ylabel("Ec / J")

% %% Potential
% figure(3)
% subplot(4, 1, 3)
h = L * (1 - cos(theta));
U = M * g * h;
plot(t, U)
title("Energia Potencial ao longo do tempo")
xlabel("t / s")
ylabel("Ep / J")
%%
% %% Mechanical
% figure(4)
% subplot(4, 1, 4)
E = K + U;
plot(t, E)
% axis([0 30 183 190])
title("Energia Mecânica ao longo do tempo")
xlabel("t / s")
ylabel("Em / J")