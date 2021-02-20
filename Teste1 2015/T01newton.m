function [XZ,n] = T01newton(funcao,funcaoprime,xa,deltac);
%
%

n=0;
delta=1;
while delta>deltac %o que traduz esta condição?
    x=xa-funcao(xa)/funcaoprime(xa);  
    n=n+1;
    delta=abs(xa-x)/abs(x);
    xa=x; %porquê?
end
XZ=x;
end