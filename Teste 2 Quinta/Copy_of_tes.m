%%%
%%%

%%%%%%%%%%%%% MUDAR QUIVER DAS MOLAS PARA PLOTS



close all
clear all
clc



%% Constantes e variaveis dadas
g=10; %m/s^2

m=1; %kg
L=1; %m
k=4; %N/m

Tf=7; %s

n=10; % numero de passos


dt=0.1;
t=0:dt:n*Tf;

%% Teto:
i = 1:length(t);
yT(i)=cos(2*pi*i*dt/Tf);




%% Euler-Cromer AZUL

% valores iniciais

yB(1) = -L;
vB(1) = Fres(g,k,yT(1),yB(1),m)*dt/m;

for i = 1:length(t)-1

    if (yB(i)>=yT(i))&& (vB(i)>0)
       vB(i)=-vB(i); 
    end
    vB(i+1) = vB(i)+Fres(g,k,yT(i),yB(i),m)*dt/m;
    yB(i+1) = yB(i)+vB(i+1)*dt;
end


vB_EC=vB;
yT_EC=yT;
yB_EC=yB;


clearvars vB
clearvars yB




Em_EC= m/2.*vB_EC.*vB_EC + m*g.*yB_EC;

%% Verlet  VERMELHO


yB(1) = -L;
vB(1) = Fres(g,k,yT(1),yB(1),m)*dt/m;


yB(2) = yB(1)+vB(1)*dt;
vB(2) = vB(1)+Fres(g,k,yT(2),yB(2),m)*dt/m;



for i = 2:length(t)-1
     

    yB(i+1) = 2*yB(i)-yB(i-1)+Fres(g,k,yT(i),yB(i),m)*dt*dt/m;
    vB(i)=(yB(i+1)-yB(i-1))/(2*dt);
    if (yB(i+1)>=yT(i-1))&& (vB(i)>0)
       vB(i)=-vB(i);
       yB(i+1)=yB(i);
       %%yB(i+1) = 2*yB(i)-yB(i-1)+Fres(g,k,yT(i),yB(i),m)*dt*dt/m;
    end
end


vB_V=vB;
yT_V=yT;
yB_V=yB;

clearvars vB
clearvars yB

vB_V(length(vB_V)+1)=vB_V(length(vB_V));

Em_V= m/2.*vB_V.*vB_V + m*g.*yB_V ;

%%%% ERRO 
maxV=max(yB_V)
maxEC=max(yB_EC)
minV=min(yB_V)
minEC=min(yB_EC)

%% Animacao

rB=0.1; %raio da bola
delta_graf=-min(yB_EC)+2;

figure(1)

for i = 1:length(t)
    cla reset
    xlim([-delta_graf/2 delta_graf/2])
    ylim([1-delta_graf 1])
    hold on
    
    % Animacao Euler Cromer (azul)
    quiver(-1,yT_EC(i),2.25,0,'b.','LineWidth',3); 
    hold on
    quiver(0,yT_EC(i),0,-(yT_EC(i)-yB_EC(i))-0.1,'b.','LineWidth',2);
    rectangle('Position',[-rB yB_EC(i)-rB 2*rB 2*rB],'Curvature',[1 1],'FaceColor','b');
    hold on
    
    % Animacao Verlet (vermelho)
    quiver(-1,yT_V(i),2.25,0,'r.','LineWidth',3); 
    hold on
    quiver(0,yT_V(i),0,-(yT_V(i)-yB_V(i))-0.1,'r.','LineWidth',2);
    rectangle('Position',[-rB yB_V(i)-rB 2*rB 2*rB],'Curvature',[1 1],'FaceColor','r'); 
    
    %pause(dt);
    drawnow;
end


%% Gráfico do yB(t)

figure(2)
i2 = 0:length(t)-1;
plot(i2*dt,yB_EC,'b-');
hold on
plot(i2*dt,yB_V,'r-');

%% Gráfico da EM

figure(3)
i2 = 0:length(t)-1;
plot(i2*dt,Em_EC,'b-');
hold on
plot(i2*dt,Em_V,'r-');







%% Forças:

%Força resultante na bola
function f=Fres(g,k,yT,yB,m)
f=-m*g+Fel(k,yT,yB);                 
end
% forca elastica da mola
function f=Fel(k,yT,yB)
f=k*abs(yT-yB);
end

