%% Constantes
close all
clear all
clc



i=0:5;
x(1) = 10;
x(2) = 20;
x(3) = 60;
x(4) = 70;
x(5) = 80;
x(6) = 81;
t = 0:0.1:5;

dt = 0.1;

interp = spline(i, x, t);
plot(t, interp, 'k-')
hold on




%% Variação da Aceleração

coef = polyfit(t, interp, 4);

x = coef(5) + coef(4) * t + coef(3) * t .^2 + coef(2) * t .^3 + coef(1) * t.^4;

v = coef(4)+ 2 * coef(3) * t + 3 * coef(2) * t .^2 + 4 * coef(1) * t.^3;

ac= 2 * coef(3)+ 6 * coef(2) * t + 12 * coef(1) * t .^2;

plot(t,x,'b-')
plot(t,v,'r-')
plot(t,ac,'g-')
legend("Spline do Moviemnto", "Posição/Tempo", "Velocidade/Tempo", "Aceleração/Tempo", 'Location', 'NorthWest')