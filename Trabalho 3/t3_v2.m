clear all
close all
clc


%% Constants
cx = 0;
cy = 0;
G = 6.67 * 10 ^-11;
c = 3 * 10 ^8;

% M - mass of the black hole (in this case mass of the M87 BH)
M = 11 * 1.98 * 10 ^30;

% Short Sealed Radius - Radius of the Event Horizon
rs = (2 * G * M) / c ^2;

%% Photons

photon = [0; 1.5 * rs];
v = [-c; 0];

%% Equations

p = 1000;
tf = 5 * rs / c;
dt = tf / p;

x = photon(1, 1);
y = photon(2, 1);

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
    
    theta = atan(photon(1, 1) / photon(2, 1));
    r = sqrt(photon(1, 1) ^2 + photon(2, 1) ^2);
    at = G * M / (r ^2) * cos(theta);
    dTheta = at * (dt / c);
    
    dTheta = abs(mod(dTheta, 2 * pi));
    theta = (mod(theta, 2 * pi))
    
    theta = theta + dTheta;
    theta = mod(theta, 2 * pi);
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    v = R * v;
    
    photon = photon + v * dt;
    
    plot(photon(1, 1), photon(2, 1), 'ob')
    drawnow limitrate nocallbacks
    axis equal
    hold off
end







