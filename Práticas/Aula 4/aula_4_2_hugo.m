close all
clear all

s0 = 'F++F++F';

for i = 1:6
    s1 = replace(s0,'F','F-F++F-F');
    s0 = s1;
end

teta = pi/3;
mrot = [cos(teta) sin(teta) ; -sin(teta) cos(teta)];
aresta = [1;0];
v=[0;0];

for i = 1:length(s1)
    if s1(i) == 'F'
        novoponto = [v(1,end);v(2,end)]+aresta;
        v = [v, novoponto];
    end
    
    if s1(i) == '+'
        aresta = mrot*aresta;
    end
    
    if s1(i)=='-'
        aresta = transpose(mrot)*aresta;
    end
end

plot(v(1,:),v(2,:),'-')
axis equal