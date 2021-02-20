clear all
close all
clc

% x = x0;
% y = 10;
% v0 = 0;
alpha = 1;
x = 0.5;
t = linspace(0, 1, 1000);
% x = linspace(-2, 2, 100);
% y = sup(alpha, x);
% plot(x, y, 'b');
% hold on


x = 0.1;
y = 0.1;
r = [x; y];
n = 1000;
nc = 0.000001;

while n > nc
   J = [];
   
    
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