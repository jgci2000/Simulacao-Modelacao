function [ md,k ] = goldensearch( ff,a,b,sw,deltac )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

k=0;
delta=1;
rate=(((sqrt(5)+1)/2)-1);

while delta>deltac
   d=rate*(b-a);
   x1=a+d;
   x2=b-d;
   if sw*ff(x2)>sw*ff(x1)
       a=x2;
       xmin=x1; 
   else
       b=x1;
       xmin=x2;
   end
   delta = abs((b - a) / xmin);
   k=k+1;
end

md=xmin;
end

