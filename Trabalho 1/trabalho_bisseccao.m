clear all
close all
clc

teta0 = 0;
teta1 = pi / 2;
deltac = 0.000000001;
tetaa = teta0;
tetab = teta1;
delta = 1;
v0 = 20;
g = 9.8;
w = 4.65;
A = 1;


while delta > deltac
    
    teta = (tetaa + tetab) / 2;
    
    if(f(tetaa,v0,w,A) * f(teta,v0,w,A) < 0)
        delta = abs(tetab - teta) / abs(teta);
        tetab = teta;
    else
        delta = abs(tetaa - teta) / abs(teta);
        tetaa = teta;
    end
end

disp(teta);
tc=10/(v0*cos(teta));

for t=0:0.01:tc
    
    x = (v0 * cos(teta)) * t;
    y = (v0 * sin(teta)) * t - 0.5 * g * t ^ 2;  
    plot(x, y, '.b', 'MarkerSize', 20);
    axis([0 20 0 12]);
    hold on
    
    p = 5 + A * cos(w * t + (pi/2));
    xf = 10;
    plot(xf, p, '.r', 'MarkerSize', 20);
    pause(0.05);
    hold off
    
end

function fx = f(ang,v0,w,A)
fx = - 5 + (10 * tan(ang) - 0.5 * 9.8 * 100 / (cos(ang) * cos(ang) * 2 * v0 ^ 2)) - A * cos(w * 10 / (v0 * cos(ang)) + pi / 2);
end