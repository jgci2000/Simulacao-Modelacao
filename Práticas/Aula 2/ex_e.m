clear all
close all

%CRIAR O POLIGNO DE n LADOS
n = 5;
ang = 2 * pi / n;
res = zeros(2, n);
p = [sqrt(2); sqrt(2)];
R = [cos(ang), -sin(ang); sin(ang), cos(ang)];

for j = 1:n+1
    p = R * p;
    res(:, j) = p;
end

%PLOT DO POLIGNO
% xlim([-40 40])
% ylim([-40 40])
% plot(res(1, :), res(2, :), '-r')

%POR O POLIGNO A RODAR EM TORNO DE UMA CIRCUNFERENCIA DE RAIO 2
%VELOCIDADE ANGULAR - W = 2PI / T
res_b1 = res;
r = 4;
ang2 = (pi / 4);
R2 = [cos(ang2), -sin(ang2); sin(ang2), cos(ang2)];
for i = 1:1000
    ang3 = i * pi / 20;
    R3 = [cos(ang3), -sin(ang3); sin(ang3), cos(ang3)];
    res = res_b1 + r;
    res = R3 * res;
    
    res = R2 * res;
    plot(res(1, :), res(2, :), '-b')
    xlim([-10 10])
    ylim([-10 10])
    pause(0.01);
    
end

%POR O POLIGNO A RODAR SOBRE SI
% res_b2 = res;
% ang2 = (pi / 4);
% R2 = [cos(ang2), -sin(ang2); sin(ang2), cos(ang2)];
%for i = 1:10
    %res = R2 * res;
    %plot(res(1, :), res(2, :), '-b')
    %pause(0.1);
%end

