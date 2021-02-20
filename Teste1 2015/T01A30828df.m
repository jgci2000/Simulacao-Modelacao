function [ dr ] = T01A30828df( t )
%T01A30828F Derivada da funcao diferenca. Teste 1 Ex. 2
%   Tiago A. Oliveira

dt=0.001;

dr=(T01A30828f(t+dt)-T01A30828f(t))/dt;
end

