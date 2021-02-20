clear all
close all

% Definition of masses

M = 1; % Massa das bolas 
Mc = 2; % Massa da caixa
L=5; % largura da caixa
%
K = 2; % constante da mola

NT = 5000; % number of time intervals
%Initial position of masses
%
Xc = 0;
%
fe= 1/3; % a value between 0 and 0.9
fd = 1/2; % a value between 0.1 and 0.9
%
Xei = Xc + fe*3*L;
Xdi  = Xei + fd*(3*L -Xei);

Xi = [Xei-L, Xdi-2*L, 1]'; %vector posicao inicial nas vari?veia transformadas;
Vi = [3, -6, 0]'; %velocidade inicial;
Vcm = 0; % velocidade inicial do centro de massa
if Vcm == 0
    Vi = Vi - (M*Vi(1,1) + M*Vi(2,1) + Mc*Vi(3,1))/(2*M+Mc);
end
%%
% Passo 2: Modos normais 

MSyst = [2*K/M, -K/M, -K/M; -K/M, 2*K/M, -K/M; -K/Mc, -K/Mc, 2*K/Mc];

[V, D] = eig(MSyst);

% Definir dt
omega_max = max(eig(MSyst));
Tmin= 2*pi/omega_max;
dt = Tmin/50;

% Passo 1: M?todo de Euler-Cromer (nas vari?veis transformaas)

V= NaN(3,NT);
X = V;
V(:,1) = Vi;
X(:,1) = Xi;

for i = 2:NT
    F = -MSyst*X(:,i-1);
    V(:,i) = V(:,i-1) +dt*F;
    X(:,i) = X(:,i-1) + dt*V(:,i);
    if(X(2,i) + L - X(1,i) < 0 || X(2,i) + 2*L > X(3,i) +3*L || X(1,i) + L <  X(3,i))
        mensagem = 'not a physical situation'
        return
    end
end
  X = X + [L; 2*L; 0]; % passar ?s vari?veis originais
   
% Energia total

ECbolas = 1/2*M*sum(V(1:2,:).^2);
ECcaixa = 1/2*Mc* V(3,:).^2;
EP = 1/2*K*( (X(1,:) - X(3,:) -L).^2 + ( X(2,:) - X(1,:) -L).^2 + (X(2,:) + X(3,:) -2*L).^2);
    
Etot = ECbolas + ECcaixa + EP;

figure
plot(0:dt:(NT-1)*dt,Etot)

figure
plot(0:dt:(NT-1)*dt,X')