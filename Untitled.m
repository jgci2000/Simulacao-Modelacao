p = [1; 2];
v = [4; 0];

n = 3;
ang = 2 * pi / n;
R = [cos(ang), -sin(ang); sin(ang), cos(ang)];
res = zeros(2, n+1);
for i = 1:n + 1
    p = R * p;
    res(:, i) = p;
end
 plot(res(1,:), res(2, :), 'ob')