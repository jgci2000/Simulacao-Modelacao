%% Constantes
clear all
close all
clc

m = 1;
L = 1;
K = 4;
g = 10;

n = 10; % número de períodos
dt = 0.1;
Tf = 7;
ti = 0;
tf = Tf * n;

t = ti:dt:tf;
i = 1:length(t);
yt(i) = cos(2 * pi * t / Tf);

%% Euler-Cromer

vb = zeros(1, length(t));
yb = zeros(1, length(t));

yb(1) = -L;
vb(1) = Fr(g, m, K, yt(1), yb(1)) * dt / m;

for i = 1:length(t) - 1
    vb(i + 1) = vb(i) + dt * Fr(g, m, K, yt(i), yb(i)) / m;
    yb(i + 1) = yb(i) + dt * vb(i + 1);
end

yb_ec = yb;
vb_ec = vb;


%% Verlet

vb = zeros(1, length(t));
yb = zeros(1, length(t));

yb(1) = -L;
vb(1) = Fr(g, m, K, yt(1), yb(1)) * dt / m;

yb(2) = yb(1) + dt * vb(1);

for i = 2:length(t) - 1
    yb(i + 1) = (2 * yb(i) - yb(i - 1) + (Fr(g, m, K, yt(i), yb(i) * dt ^2) / m)) * (dt ^2);
    vb(i) = (yb(i + 1) - yb(i - 1)) / (2 * dt);
    
    if (yb(i+1)>=yt(i-1))&& (vb(i)>0)
       vb(i)=-vb(i);
       yb(i+1)=yb(i);
       %%yB(i+1) = 2*yB(i)-yB(i-1)+Fres(g,k,yT(i),yB(i),m)*dt*dt/m;
    end
end

yb_v = yb;
vb_v = vb;

%% Plots

figure(1)
plot(t, yb_ec, '-b')
hold on
plot(t, yb_v, '-r')
legend("Euler-Cromer y", "Verlet y")

figure(2)
plot(t, vb_ec, '-b')
hold on
plot(t, vb_v, '-r')
legend("Euler-Cromer v", "Verlet v")

%% Energia Mecânica

Em_ec = 0.5 * m * vb_ec .^2 + m * g * yb_ec;

Em_v = 0.5 * m * vb_v .^2 + m * g * yb_v;

figure(3)
plot(t, Em_ec, '-b')
hold on
plot(t, Em_v, '-r')
legend("Euler-Cromer Em", "Verlet Em")

%% Forças

function f1 = Fr(g, m, K, yt, yb)
f1 = - m * g + Fel(K, yt, yb);
end

function f2 = Fel(K, yt, yb)
f2 = K * abs(yt - yb);
end