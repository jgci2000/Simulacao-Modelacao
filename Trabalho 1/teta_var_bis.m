close all
clear all
clc

for v0 = 13:200
    teta0 = 0;
    teta1 = pi/2;
    deltac = 0.000000001;
    tetaa = teta0;
    tetab = teta1;
    delta = 1;
    g = 9.8;
    w = 4.65;
    A = 1;

    while delta > deltac
        teta = (tetab + tetaa) / 2;
        if(f(tetaa, v0, w, A) * f(teta, v0, w, A) < 0)
            delta = abs(tetab - teta) / abs(teta);
            tetab = teta;
        else
            delta = abs(tetaa - teta) / abs(teta);
            tetaa = teta;
        end
    end
    
    plot(v0, teta, '.r');
    axis([0 200 0 2]);
    title("teta em função de v0 - Método da Bissecção");
    xlabel("v0");
    ylabel("teta");
    hold on

end

function fx = f(ang,v0,w,A)
fx = - 5 + (10 * tan(ang) - 0.5 * 9.8 * 100 / (cos(ang) * cos(ang) * 2 * v0 ^ 2)) - A * cos(w * 10 / (v0 * cos(ang)) + pi / 2);
end