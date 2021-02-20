clear all
close all
clc

axis ([-10 30 0 30]);
hold on

% parametros alterar
mpend = 4 % massa do pendulo
mbala = 0.1 % massa da bala
vinicial = 300  % velocidade inicial da bala
 
%pendulo
xp = 10;
yp = 5;

% teto
yteto = 18;
tetox = [5, 15];
tetoy = [yteto, yteto];

L = yteto-yp; %comprimento do fio

%bala 
xb = -5;
yb = yp;

g = 9.8;
t= 0;
dt = 0.05;

theta = 0;


% movimento 

while xb < xp
    plot(tetox, tetoy, 'k', 'Linewidth',10); %plot teto
    hold on
    plot([xp xp],[yp yteto], 'k', 'Linewidth', 2); %plot fio
    plot(xp, yp, '.r','markersize', 120); %plot pendulo
    plot(xb, yb,'.','markersize', 20); %plot bala
    hold off
    yb=yp;
    xb= + vinicial*t;
    axis ([-10 30 0 30]);
    pause(0.005);
    t = t+0.005;
end

yb = yp;
xb = xp;

v = (mbala/(mbala+mpend))* vinicial %velocidade após a colisão (conservação do momento linear)
theta_f = acos(1-((v*v)/(2*g*L))) %angulo maximo (conservação da enegia mecanica)
vang = sqrt(g/L);
hmax = (v.^2)/(2*g)
   
for t = t:dt:50
    theta = theta_f*cos(vang*t-pi/2);
    x = xp + L*sin(theta);
    y = yp + L*(1-cos(theta));
    plot(tetox, tetoy, 'k', 'Linewidth',10); %plot teto
    hold on
    plot([x xp],[y yteto], 'k', 'Linewidth', 2); %plot fio
    plot(x, y, '.r','markersize', 120); %plot pendulo
    plot(x, y,'.','markersize', 20); %plot bala
    hold off

    axis ([-10 30 0 30]);
    pause(dt);
end