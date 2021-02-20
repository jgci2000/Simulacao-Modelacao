clear all
close all

n = 7;
ang = (2 * pi) / n;
R = [cos(ang), -sin(ang); sin(ang), cos(ang)];
p = [2; 2];
v = [sqrt(2); sqrt(2)];
res = zeros(2, n);
res(:, 1) = p;

for i = 1:n + 1
    p = v + p;
    v = R * v;
    res(:, i) = v;
end

plot(res(1, :), res(2, :), '-r')
axis equal