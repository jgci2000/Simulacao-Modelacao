function [ y ] = tffiting( T )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global x t

xx=5+cos((2*pi./T)*t);  % Funcao a 'fitar'

y=sum((xx-x).^2);   % Minimo do residual e a solucao
end

