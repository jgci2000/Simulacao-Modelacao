close all
clear all

% CRIAR O POLIGNO DE n LADOS
n = 4;
ang = 2 * pi / n;
res = zeros(2, n);
p = [3; 3];

R = [cos(ang), -sin(ang); sin(ang), cos(ang)];

for j = 1:n+1
    p = R * p;
    res(:, j) = p;
end

% PLOT DO POLIGNO
plot(res(1, :), res(2, :), '-r', 'LineWidth', 1.25)
axis equal
hold on

% TRANSLA플O

res_t = res;
res_t = res_t + 10;
plot(res_t(1, :), res_t(2, :), '-b')

%ROTA플O

res_r = res;
ang_r = pi / 4;
R_r = [cos(ang_r), -sin(ang_r); sin(ang_r), cos(ang_r)];
res_r = R_r * res_r;
plot(res_r(1, :), res_r(2, :), '-k')

% AMPLIA플O/REDU플O

res_a = res;
res_red = res;


% DISTOR플O

res_d = res;

