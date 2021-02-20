clear all
close all
clc


A = 1;
w = 4.65;
v0 = 20;
g = 9.8;

deltac = 0.000000001;
teta0 = 0;
teta1 = pi / 2;
tetaa = teta0;
tetab = teta1;
delta = 1;

lambda = 0.1;
f1 = f(tetaa, v0, g, A, w);
teta = tetaa - lambda .* f1;
f2 = f(teta, v0, g, A, w);

if abs(f1) < abs(f2)
    lambda = -lambda;
end

while delta > deltac
    f1 = f(teta, v0, g, A, w);
    teta = teta - lambda .* f1;
    delta = abs(lambda .* f1) / abs(teta);
end

disp(teta);
tc = 10 / (cos(teta) * v0);

for t = 0:0.01:tc
    
    x = (v0 * cos(teta)) * t;
    y = (v0 * sin(teta)) * t - 0.5 * g * t ^ 2;
    plot(x, y, '.b', 'MarkerSize', 20);
    axis([0 20 0 12]);
    hold on
    
    p = 5 + cos(w * t + (pi / 2));
    xf = 10;
    plot(xf, p, '.r', 'MarkerSize', 20);
    drawnow update
%     pause(0.05);
    hold off
    
end



function f1 = f(ang, v0, g, A, w)
f1 = - 5 + (10 * tan(ang) - 0.5 * g * (100 / (cos(ang) * cos(ang) * v0 ^ 2))) - A * cos(w * 10 / (v0 * cos(ang)) + pi / 2);
end