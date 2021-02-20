close all
clear all

%CRIAR O POLIGNO DE n LADOS
n = 5;
ang = 2 * pi / n;
res = zeros(2, n);
p = [3; 1];

R = [cos(ang), -sin(ang); sin(ang), cos(ang)];

for j = 1:n+1
    p = R * p;
    res(:, j) = p;
end

%PLOT DO POLIGNO
plot(res(1, :), res(2, :), '-r')
axis ([0 50 -1 50])
hold on

%ROTAÇÃO DE 45 E MOVIEMENTO DE 5 UNIDADES, DEPOIS RODA 30 GRAUS
res_b = res;
ang3 = pi / 6;
R3 = [cos(ang3), -sin(ang3); sin(ang3), cos(ang3)];
for i = 1:5
    ang2 = (i * pi / 4);
    R2 = [cos(ang2), -sin(ang2); sin(ang2), cos(ang2)];
    res = R2 * res_b;
    res = res + (5 * i);
    res = R3 * res;
    plot(res(1, :), res(2, :), '-b');
end


