clear all
close all
clc
%A cada F tem de desenhar 5 unidades
%A cada + tem de rodar +60º
%A cada - tem de rodar -60º

s0 = 'F++F++F';

for i = 1:7
    s1 = replace(s0, 'F', 'F-F++F-F');
    s0 = s1;
end

ang = pi / 3;
R_p = [cos(ang), -sin(ang); sin(ang), cos(ang)];
R_n = [cos(-ang), -sin(-ang); sin(-ang), cos(-ang)];
aresta = [5; 0];
s1_length = length(s1);
v = [0; 0];

for i = 1:s1_length
    if s1(i) == 'F'
        np = [v(1, end); v(2, end)] + aresta;
        v = [v, np];
    end
    
    if s1(i) == '-'
        aresta = R_n * aresta;
    end
    
    if s1(i) == '+'
        aresta = R_p * aresta;
    end
end

plot(v(1,:), v(2,:), '-b');
axis equal