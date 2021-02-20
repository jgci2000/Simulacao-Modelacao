clear all
close all
clc


figure(1)

x = linspace(0, 20, 21);
lambda = 7;

plot(x, f(x, lambda), 'or');
hold on

xpol = linspace(-10, 30, 100);
ypol = interp1(x, f(x, lambda), xpol);
ypol1 = interp1(x, f(x, lambda), xpol, 'PCHIP');
ypol2 = interp1(x, f(x, lambda), xpol, 'SPLINE');

plot(xpol, ypol, '.b');
plot(xpol, ypol1, 'ok');
plot(xpol, ypol2, 'oy');

legend('f(x)', 'ypol', 'ypol1', 'ypol2', 'Location', 'northwest');


figure(2)

plot(x, f(x, lambda), 'or');
hold on

a = polyfit(x, f(x, lambda), 2);
y1 = a(3) + a(2) * x + a(1) * x .^ 2;
plot(x, y1, 'ok')

legend('f(x)', 'polynomial fit', 'Location', 'northwest');


function y = f(x, lambda)
y = exp(x / lambda);
end