%% Polar Coordinates - Trabalho 3
clear all
close all
clc

%% Constants
cx = 0;
cy = 0;
G = 6.67 * 10 ^-11;
c = 3 * 10 ^8;

% M - mass of the black hole
% 11 times the mass of the sun
M = 11 * 1.98 * 10 ^30;

% Short Sealed Radius - Radius of the Event Horizon
rs = (2 * G * M) / c ^2;

%% Photons
% (r, theta)

photon = [4 * rs; pi / 3];

%% Equations

p = 100;
tf = 2 * rs / c;
dt = tf / p;

v = c;

r = photon(1, 1);
theta = photon(2, 1);

for t = linspace(0, tf, p)
    %% Black Hole
    
    % Acretion Disk
    Circle(cx, cy, rs * 3, [0.5 0.5 0.5], 20)
    hold on
    % Photon Sphere
    Circle(cx, cy, 1.5 * rs, [1 0.40 0.05], 10)
    % Center
    FilledCircle(cx, cy, rs, 'k')
    
    %% Photon Movement
    
    at = G * M / (r ^2) * cos(theta);
    
    dTheta = at * (dt / c);
    dTheta = mod(dTheta, 2 * pi);
    
    theta = theta + dTheta;
    
    if theta > 2 * pi 
        theta = mod(dTheta, 2 * pi);
    end
    
    dr = (v + r * dTheta) * dt;
    r = r + dr;
    
    x = r * cos(theta);
    y = r * sin(theta);
    
    
    plot(x, y, 'ok');
    axis equal
%     axis([-2 * 10 ^5, 3 * 10 ^5, -2 * 10 ^5, 2 * 10 ^5])
%     pause
    drawnow update
    hold off
end







