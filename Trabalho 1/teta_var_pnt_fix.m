clear all
close all
clc


for v0 = 13:200
    
    A = 1;
    w = 4.65;
    g = 9.8;

    deltac = 0.000000001;
    teta0 = 0;
    teta1 = pi / 2;
    tetaa = teta0;
    tetab = teta1;
    delta = 1;

    lambda = 0.1;
    f1 = f(tetaa, v0, g, A, w);
    teta = tetaa - lambda .* f1;
    f2 = f(teta, v0, g, A, w);

    if abs(f1) < abs(f2)
        lambda = -lambda;
    end

    while delta > deltac
        f1 = f(teta, v0, g, A, w);
        teta = teta - lambda .* f1;
        delta = abs(lambda .* f1) / abs(teta);
    end
    
    plot(v0, teta, '.r');
    axis([0 200 0 2]);
    title("teta em função de v0 - Método do Ponto Fixo");
    xlabel("v0");
    ylabel("teta");
    hold on
end


function f1 = f(ang, v0, g, A, w)
f1 = - 5 + (10 * tan(ang) - 0.5 * g * (100 / (cos(ang) * cos(ang) * v0 ^ 2))) - A * cos(w * 10 / (v0 * cos(ang)) + pi / 2);
end