%%%
%%%

close all
clear all
clc

%%%% t=0 ----> i=1 .....

i=0:5;
x(1)=10;
x(2)=20;
x(3)=60;
x(4)=70;
x(5)=80;
x(6)=81;
t=0:0.1:5;

interp=spline(i,x,t);
plot(t,interp,'k-')
hold on




%% variacao da Aceleracao nao constante
a=polyfit(t,interp,4);

eq_x= a(5)+ a(4)*t+ a(3)*t.*t+ a(2)*t.^3 +a(1)*t.^4 ;

eq_v= a(4)+ 2*a(3)*t+ 3*a(2)*t.*t + 4*a(1)*t.^3 ;

eq_ac= 2*a(3)+ 6*a(2)*t + 12*a(1)*t.^2 ;

plot(t,eq_x,'b-')
hold on
plot(t,eq_v,'r-')
hold on
plot(t,eq_ac,'g.')
hold on

%%%%Aumento da precisao
dt2=0.001;
t2=0:dt2:5;
ac2=ac1(a,t2);
v2=v1(a,t2);
x2=x1(a,t2);

%%%% a) Acelaraco minima
iteracao=find(ac2==min(ac2));
t_ac_min=(iteracao*dt2)-dt2 % instante da aceleracao minima
x_ac_min=x2(iteracao) %x onde acelaercao minima
ac_min=ac2(iteracao) %aceleracao minima
%%%ac_min2=min(ac2


%%%% b) imobilizado
c=(1/dt2)/100;

iteracaoI=find(round(v2*c)==0);
t_imobilizado=(iteracaoI*dt2)-dt2 % instante para

x_imobilizado=x2(iteracaoI) %x onde para






function f=x1(a,t)
f= a(5)+ a(4)*t+ a(3)*t.*t+ a(2)*t.^3 +a(1)*t.^4 ;
end
function f=v1(a,t)
f= a(4)+ 2*a(3)*t+ 3*a(2)*t.*t + 4*a(1)*t.^3 ;
end
function f=ac1(a,t)
f= 2*a(3)+ 6*a(2)*t + 12*a(1)*t.^2 ;
end



%%%%%%%%%%%%%%%%%%%%--------------%%---------------%%%%%%%%%%%%%%%%%%%%%%%%

% % % % %% Aceleracao constante
% % % % a=polyfit(t,interp,2);
% % % %
% % % % eq_x= a(3)+ a(2)*t+ a(1)*t.*t;
% % % %
% % % % eq_v= a(2)+ 2*a(1)*t;
% % % %
% % % % eq_ac= 2*a(1);
% % % %
% % % % plot(t,eq_x,'b-')
% % % % hold on
% % % % plot(t,eq_v,'r-')
% % % % hold on
% % % % plot(t,eq_ac,'g.')
% % % % hold on
% % %
% % % %% Aceleracao nao constante
% % % a=polyfit(t,interp,3);
% % %
% % % eq_x= a(4)+ a(3)*t+ a(2)*t.*t+ a(1)*t.^3;
% % %
% % % eq_v= a(3)+ 2*a(2)*t+ 3*a(1)*t.*t ;
% % %
% % % eq_ac= 2*a(2)+ 6*a(1)*t;
% % %
% % % plot(t,eq_x,'b-')
% % % hold on
% % % plot(t,eq_v,'r-')
% % % hold on
% % % plot(t,eq_ac,'g.')
% % % hold on