clear all;
close all;

axis ([0 1 0 1]);

N = 3000;
x = rand(1,N); %coordenada aleatórias
y = rand(1,N);
dmin=100;
flag = zeros(1,N); %pontos ligados
dt = 0.05;

plot(x,y,'.k','markersize',10); %plot pontos
hold on

%calcula a menor distancia e une os pontos com menor distância entre eles
for i=1:N
    for j=i+1:N
        dist = sqrt(((x(i)-x(j)).^2)+ ((y(i)-y(j)).^2));
        if dmin>dist
            dmin=dist;
            primeiro_ponto=i;
            segundo_ponto=j;
        end
    end
end
flag(1, primeiro_ponto) = 1;
flag(1, segundo_ponto) = 1;
xx = [x(1,primeiro_ponto), x(1,segundo_ponto)];
yy = [y(1,primeiro_ponto), y(1,segundo_ponto)];
plot(xx,yy,'k','Linewidth',1);
hold on
pause(dt);

%arvore de menor extensão - animação
for L=2:N-1
    dmin=100;
    for i=1:N
        for j=i+1:N
            if flag(1,i) ~= flag(1,j)
                dist = sqrt(((x(i)-x(j)).^2)+ ((y(i)-y(j)).^2));
                if dmin>dist
                    dmin=dist;
                    primeiro_ponto=i;
                    segundo_ponto=j;
                end
            end
        end
    end
    flag(1, primeiro_ponto) = 1;
    flag(1, segundo_ponto) = 1;
    xx = [x(1,primeiro_ponto), x(1,segundo_ponto)];
    yy = [y(1,primeiro_ponto), y(1,segundo_ponto)];
    plot(xx,yy,'k','Linewidth',1);
    pause(dt);
    %drawnow
end
            



