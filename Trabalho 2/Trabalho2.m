clear all
close all
clf

%%
% 1º passo
mov = VideoReader('bola.mp4');
video=read(mov);
s = size(video);
nFrames = s(4);
framerate = mov.FrameRate;


%%
% 2º passo

%mudança de cor
frame1 = video(:,:,1,1);
colormap gray
frame1(frame1<=45)=1;
frame1(frame1~=1)=0;
imagesc(frame1);


%3º passo
[xp,yp]=find(frame1==1);
centro = mean([xp,yp]);
hold on
axis equal
ylabel("Altura")
title("Primeiro Frame c/ Centro de Massa")
plot(centro(2),centro(1),'or')

%%
%4ºpasso
posY=1-zeros(1,866);
frames = zeros(1,866);

for i=1:nFrames
    
    frames(i)=i;
    frame2=video(:,:,1,i);
    frame2(frame2<=45)=1;
    frame2(frame2~=1)=0;
    imagesc(frame2);
    [xp,yp]=find(frame2==1);
    centro = mean([xp,yp]);
    posY(1,i)=680-centro(1);
    hold on
    axis equal
    plot(centro(2),centro(1),'or')
    drawnow
    hold off
    
end

%%
%5º e 6º passo

hold off
frames2s =frames./framerate;
dim=172.2310-123.9238;
pixel2m=1*dim/0.05;
posYmetros=posY./pixel2m;
plot(frames2s,posYmetros,'-b');
ylabel("Altura / m")
xlabel("Tempo / s")
title("Movimento da Bola ao Longo do Tempo")
hold on
hold off


%%
%7º passo
i=1:866;
dt=1/framerate;
t=i.*dt;
tempos = linspace(0,6,6*framerate);
interpolacao = spline(frames2s,posYmetros,tempos);

plot(tempos,interpolacao)
xlabel("Tempo / s")
ylabel("Altura / m")
title("Spline")

[maximos,t_max] = findpeaks(interpolacao,tempos);
[minimos,t_min] = findpeaks(-interpolacao,tempos);

%%
%8º passo

yextrap = interp1(frames2s,posYmetros,linspace(0,frames2s(length(frames2s))),'spline');
plot(linspace(0,frames2s(length(frames2s))),yextrap)

ylabel("Altura / m")
xlabel("Tempo / s")
title("Interpolação")



% tg = frames2s((find(frames2s==tmin(1))):(find(frames2s==tmin(2))));
% yg = -posYmetros((find(frames2s==tmin(1))):(find(frames2s==tmin(2))));
% a=polyfit(tg,-yg,2);
% yy=polyval(a,tg);
% plot(tg,yy)
% g1=-2*a(1);
%%
for c=1:17
    g2(c)=8*maximos(c)/((t_max(c)-t_min(c))^2);
end
    mean(g2)
%%
%parte2 - metodo de Euler

g=-9.8;
e=sqrt(maximos(2)/maximos(1));
dt=0.0001;
v(1)=0;
y(1)=1;
t= 0:dt:6;

for i = 1:length(t)-1
    if y(i)<=0 && v(i)<0
        v(i+1)=v(i)*(-e);
    else
        v(i+1) = v(i)+g*dt;
    end
    y(i+1) = y(i)+v(i)*dt;
end

figure(1)
plot(t,y,'r-')
title("Método de Euler")
xlabel("Tempo / s")
ylabel("Altura / m")

%%   
%parte2 - metodo de  Euler-Cromer

g=-9.8;
e=sqrt(maximos(2)/maximos(1));
dt=0.0001;
v(1)=0;
y(1)=1;
t= 0:dt:6;

for i = 1:length(t)-1
    if y(i)<=0 && v(i)<0
        v(i+1)=v(i)*(-e);
    else
        v(i+1) = v(i)+g*dt;
    end
    y(i+1) = y(i)+v(i+1)*dt;
end

figure(2)
plot(t,y,'b--')
title("Método de Euler-Cromer")
xlabel("Tempo / s")
ylabel("Altura / m")
