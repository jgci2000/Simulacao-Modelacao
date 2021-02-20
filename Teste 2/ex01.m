%% Constantes

close all
clear all
clc

m = 1;
L = 1;
K = 2; % O user não sabe
g = 10;

n = 10; % Número de períodos
Tf = 7;
ti = 0;
dt = 0.1;
tf = Tf * n;

t = ti:dt:tf;
F = cos(2 * pi * t / Tf) / 5;

%% Euler-Cromer

vb = zeros(1, length(t));
xb = zeros(1, length(t));

xb(1) = L;
vb(1) = Fr(K, xb(1), F(1)) * dt;

for i = 1:length(t) - 1
%     if (xb(i) <= 0 && vb(i) < 0) || (xb(i) >= 2 && vb(i) < 0)
%         vb(i) = -vb(i);
%     end
    
    vb(i + 1) = vb(i) + dt * Fr(K, xb(i), F(i)) / m;
    xb(i + 1) = xb(i) + dt * vb(i + 1);
end

xb_ec = xb;
vb_ec = vb;

%% Verlet

vb = zeros(1, length(t));
xb = zeros(1, length(t));

xb(1) = L;
vb(1) = Fr(K, xb(1), F(1)) * dt;

xb(2) = xb(1) + vb(1) * dt;

for i = 2:length(t) - 1
    xb(i + 1) = 2 * xb(i) - xb(i - 1) + dt ^2 * Fr(K, xb(i), F(i)) / m;
    vb(i) = (xb(i + 1) - xb(i - 1)) / (2 * dt);
    
%     if (xb(i + 1) >= 2 && vb(i) > 0) || (xb(i + 1) <= 0 && vb(i) < 0)
%         vb(i) = - vb(i);
%         xb(i + 1) = xb(i);
%     end
end

xb_v = xb;
vb_v = vb;

%% Plots

figure(1)
plot(t, xb_ec)
hold on
plot(t, xb_v)
legend("Euler-Cromer x", "Verlet x")

figure(2)
plot(t, vb_ec)
hold on
plot(t, vb_v)
legend("Euler-Cromer v", "Verlet v")

%% Energia Mecânica

em_ec = 0.5 * m * vb_ec .^2;
em_v = 0.5 * m * vb_v .^2;

figure(3)
plot(t, em_ec)
hold on
plot(t, em_v)
legend("Euler-Cromer Em", "Verlet Em")

% %% Forças
% 
% function f = Fr(K, x, F)
% f = 2 * Fe(K, x) + F;
% end
% 
% function f = Fe(K, x)
% f = - K * x;
% end