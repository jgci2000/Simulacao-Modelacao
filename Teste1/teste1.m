clear all
close all


t0 = 0;
t1 = 3;
ta = t0;
tb = t1;
deltac = 0.000001;
delta = 1;

if(dis(t0)*dis(t1) > 0)
    disp('Escolha outros pontos. Não há zeros entre estes dois valores');
    return;
end

while delta > deltac
    
    t = (ta + tb) / 2;
    
    if(dis(ta) * dis(t) < 0)
        delta = abs(tb - t) / abs(t);
        tb = t;
    else
        delta = abs(ta - t) / abs(t);
        ta = t;
    end
end

a = t;
disp(t);

for t = linspace(0, a + 0.01, 1000)
    
    if t == a + 0.01
        pause(5);
        teste1 %ATENÇÃO: aqui tem de estar o nome do ficheiro, caso contrárip, não funciona.
    end
    hold off
    
    x = 1;
    ye = 10 + 10 * t - 0.5 * 9.8 * t ^ 2;
    yo = 2 + cos(t);
    
    plot(x, ye, '.r', 'MarkerSize', 50);
    hold on
    plot(x, yo, '.b', 'MarkerSize', 50);
    axis([0.9 1.1 0.5 16.5]);
    pause(0.001);
    
end

function f = dis(t)
f = 2 + cos(t) - (10 + 10 * t - 0.5 * 9.81 * t ^ 2);
end