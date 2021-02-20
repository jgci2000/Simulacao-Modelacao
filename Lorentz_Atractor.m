clear all
close all
clc

x = 0.01;
y = 0.01;
z = 0.01;
dt = 0.01;

a = 10;
b = 28;
c = 8.0/3.0;

pnts = zeros(3, 1201);
pnts(1, 1) = x;pnts(2, 1) = y;pnts(3, 1) = z;

for i = 2:1201
    
    dx = (a * (y - x)) * dt;
    dy = (x * (b - z) - y) * dt;
    dz = (x * y - c * z) * dt;
    
    x = x + dx; y = y + dy; z = z + dz;
    
    pnts(1, i) = x; pnts(2, i) = y; pnts(3, i) = z;
end

direction = [1 1 0];

for i = 1:1201
    p = plot3(pnts(1, i), pnts(2, i), pnts(3, i), '.r');
    hold on
    pause(0.0001);
end
