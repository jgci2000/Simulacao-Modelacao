clear all
close all

p_i = [1; 1];
res = zeros(2, 5);
res(:, 1)= p_i;

for n = 1:1:5
    
    %n do R tem de ser o número de lados do poligno...
    R = [cos(2*pi/5), -sin(2*pi/5); sin(2*pi/5), cos(2*pi/5)];
    p_f = R * p_i;
    
    res(:, n + 1) = p_f;
    p_i = p_f;
end

plot(res(1, :), res(2, :), '-r')
