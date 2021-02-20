clear all
close all
clc

x0 = 0;
y0 = 10;
v0x = 2.5;
v0y = 0;


alpha = 1;
x = 0;

lambda = 0.1;
ta = 0;
d1 = dif(alpha, x, ta);
t = ta - lambda .* d1;
d2 = dif(alpha, x, t);
delta = 1;
deltac = 0.000001;

if abs(d1) < abs(d2)
    lambda = -lambda;
end

while delta > deltac
    d1 = dif(alpha, x, t);
    t = t - lambda * d1;
    delta = abs(lambda * d1) / abs(t);
end
disp(t);

for t = linspace(0, t, 100)
    hold off
    
    x = linspace(-2, 2, 250);
    xb = x0 + v0x*t;
    yb = bola(t);
    
    plot(xb, yb, '.r', 'MarkerSize', 15);
    axis([-2 2 5 20]);
    hold on
    
    s = sup(alpha, x);
    plot(x, s);
    pause(0.01);
end



function f = dif(a, x, t)
f = - bola(t) + sup(a, x);
end

function y = bola(t)
y = 10 - 9.81 * t .^ 2;
end

function f = sup(a, x)
f = -a*x .^ 2 + x .^ 4 + 6;
end