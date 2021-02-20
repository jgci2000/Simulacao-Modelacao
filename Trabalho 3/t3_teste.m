clear all
close all
clc


%% Constants
cx = -10 ^13;
cy = 20;
G = 6.67 * 10 ^-11;
c = 3 * 10 ^8;

% M - mass of the black hole (in this case mass of the M87 BH)
M = 2.6 * 10 ^9 * 1.98 * 10 ^30;

% Short Sealed Radius - Radius of the Event Horizon
rs = (2 * G * M) / c ^2;

% Acretion Disk  - Radius = 3 * rs / Matter orbits here
% Photon Sphere - Radius = 1.5 * rs / Light orbits here

%% Photons
% Photons above 2.6 * rs, starting at blackhole.mid

n = 1;
photons = zeros(2, n + 1);
x0 = -5.5 * 10 ^13;
y0 = 0;
y = y0;
yf = 2.6 * rs;
i = 2;
photons(1, 1) = x0;
photons(1, 2) = y0;
y_inc = yf / n;

theta_photon = zeros(1, n + 1);
theta_photon(1, :) = 0;
veli = zeros(2, n + 1);
veli(1, :) = -c;


photons(1,:)=x0;
photons(2,:)=y0:y_inc:yf;

% tic
% while y <= yf
%     y = y + y_inc;
%     photons(1, i) = x0;
%     photons(2, i) = y;
%     i = i + 1;
% end
% toc
%% Plots Constants

tf = 2 * x0 / c;
p = 1000;
dt = tf / p;
% x0 = photons(1, :);
% y0 = photons(2, :);
s = size(photons);
size_photons = s(1, 2);
% r_array = zeros(1, size_photons);

force = zeros(2, size_photons);
r = zeros(1, size_photons);
at = zeros(1, size_photons);
delta_theta = zeros(1, size_photons);
theta = zeros(1, size_photons);
vel = veli;
%% Plots
for t = linspace(0, tf, p)
    %% Black Hole
    
    % Acretion Disk
    Circle(cx, cy, rs * 3, [0.5 0.5 0.5], 35)
    hold on
    % Photon Sphere
    Circle(cx, cy, 1.5 * rs, [1 0.40 0.05], 20)
    % Center
    FilledCircle(cx, cy, rs, 'k')
    
    %% Forces

    for i = 1:size_photons
        force(:, i) = photons(:, i) - [cx; cy];
%         theta(1, i) = atan(force(2, i) / force(1, i));
        r(1, i) = sqrt(force(1, i) ^2 + force(2, i) ^2);
%         theta(1, i) = atan(photons(2, i) / photons(1, i));
        theta(1, i) = asin(photons(2, i) / r(1, i))
        at(1, i) = G * M / (r(1, i) ^2);
        delta_theta(1, i) = at(1, i) * (dt / c) * sin(theta_photon(1, i) - theta(1, i));
%         delta_theta(1, i) = delta_theta(1, i) / abs(1 - ((2 * G * M) / r(1, i) * c ^2));
        theta_photon(1, i) = theta_photon(1, i) + delta_theta(1, i);
        R = [cos(theta_photon(1, i)), -sin(theta_photon(1, i)); sin(theta_photon(1, i)), cos(theta_photon(1, i))];
        vel(:, i) = R * vel(:, i);
    end
    
    %% Photons

    photons = photons + vel * dt;
    plot(photons(1, :), photons(2, :), '.y', 'Color', [1 0.05 0.05]);
    
    %% Axis
%     set(gca,'Color',[0 0 0])
        axis equal
%         axis normal
    xi = -6.5 * 10 ^13;
    xf = 3 * 10 ^13;
    yi = -5.5 * 10 ^13;
    yf = 5.5 * 10 ^13;
    axis([xi, xf, yi, yf])
    drawnow update
    hold off
end
%%







