function [ r ] = T01A30828f( t )
%T01A30828F Funcao diferenca. Teste 1 Ex. 2
%   Tiago A. Oliveira

global x01 x02 A0 b m w phi1 phi2

r1=x01+A0.*exp(-(b.*t)./(2.*m)).*sin(w.*t+phi1);
r2=x02+A0.*exp(-(b.*t)./(2.*m)).*sin(w.*t+phi2);

r=r2-r1;
end

