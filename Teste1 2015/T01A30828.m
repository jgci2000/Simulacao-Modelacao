% Simulacao e Modelacao
% (2014/2015)
% Proposta de Resolucao 1º Teste - Versao A

clear all
close all

% Declaracao das variaveis
global x01 x02 A0 b m w phi1 phi2

x01=0;
x02=1;
A0=2;
b=2;
m=2;
w=0.1;
phi1=0;
phi2=pi;

% Determinacao do momento de impacto (Met. Newton)
prec=1e-4; % Precisao

[XZ,n] = T01newton(@T01A30828f,@T01A30828df,2,prec);
% XZ converge para um valor negativo, logo os osciladores nao colidem.

% Animacao
tf=5;
t=linspace(0,tf,200);

for i=1:length(t)
    x1=x01+A0*exp(-(b*t(i))/(2*m))*sin(w*t(i)+phi1);
    x2=x02+A0*exp(-(b*t(i))/(2*m))*sin(w*t(i)+phi2);
    
    plot(x1,0.03,'bs',x2,0.03,'rs',[-3 3],[0 0],'k-')
    axis([-3 3 -1 1])
    axis off
    pause(0.01)
end