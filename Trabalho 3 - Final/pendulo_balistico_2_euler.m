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
m1 = 0.300;
v1 = 250;
xb = -10;
yb = -L;

%Constants
g = 9.81;
theta0 = acos(-((v1 / (1 + (m2 / m1))) ^2 / (2 * g * L)) + 1);

b = 0.15; % Damping Strength
M = m1 + m2;

%Plot Constants
nPoints = 1001;
dt = 0.04;
omega = zeros(nPoints, 1);
theta = zeros(nPoints, 1);
t = zeros(nPoints, 1);
theta(1) = theta0;

%% 
figure(1)
while xb < 0
    %% Bullet Animation
    tf = (x - xb) / v1;
    dt2 = 0.001;
    
    for t1 = 0:dt2:tf
        xb = xb + v1 * dt2;
        plot(xb, yb, '.r', 'MarkerSize', 20)
        hold on
        plot(x, y, '.r', 'MarkerSize', 100)
        line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
        axis([-10 10 -L-2 0])
        drawnow update
        hold off
    end
end

%% From (0, 0) to initial position

for dTheta = 0:0.05:theta0
    x = L * sin(dTheta);
    y = - L * cos(dTheta);
    
    plot(x, y, '.r', 'MarkerSize', 100)
    hold on
    line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
    axis([-10 10 -L-2 0])
    drawnow update
    hold off
end
%% ODE Solver

if b == 0
    for i = 1:nPoints - 1
        omega(i + 1) = omega(i) - (g / L) * theta(i) * dt; 
        theta(i + 1) = theta(i) + omega(i + 1) * dt;
        
        t(i + 1) = t(i) + dt;
    end
else
    for i = 1:nPoints - 1
        omega(i + 1) = omega(i) - (g / L) * theta(i) * dt - b * omega(i) * dt; % Last Expression is for the Damping Strength
        theta(i + 1) = theta(i) + omega(i + 1) * dt;
        
        t(i + 1) = t(i) + dt;
    end
end
% 1001 points equals to 40 secunds
% plot(t, theta, 'r')

%% Animation Plots

for i = 1:nPoints
    x = L * sin(theta(i));
    y = - L * cos(theta(i));
    
    plot(x, y, '.r', 'MarkerSize', 100)
    line([0 x], [0 y], 'LineWidth', 1.5, 'Color', 'red')
    axis([-10 10 -L-2 0])
    drawnow update
    hold off
end


%% Theta var funtion of time
% figure(2)
% subplot(4, 1, 1)
plot(t, theta)
title("Variação de theta ao longo do tempo")
xlabel("t / s")
ylabel("theta / rad")
%%
% %% Kinetic Energy function of time
% figure(3)
% subplot(4, 1, 2)

K = 0.5 * M * (omega .* L) .^2;
% plot(t, K)
title("Energia Cinética ao longo do tempo")
xlabel("t / s")
ylabel("Ec / J")
% %% Potential Energy function of time
% figure(4)
% subplot(4, 1, 3)

h = -(L * cos(theta) - L);
U = g * M * h;

% plot(t, U)
title("Energia Potencial ao longo do tempo")
xlabel("t / s")
ylabel("Ep / J")
% %% Mechanical Energy function od time
% figure(5)
% subplot(4, 1, 4)

E = K + U;

plot(t, E)
% axis([0 30 180 193])
title("Energia Mecânica ao longo do tempo")
xlabel("t / s")
ylabel("Em / J")

