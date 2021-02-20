function loren3
clear;clf
global A B R 
A = 10; 
B = 8/3; 
R = 28; 
u0 = 100*(rand(3,1) - 0.5); 
[t,u] = ode45(@lor2,[0,100],u0); 
N = find(t>10); v = u(N,:);

x = v(:,1);
y = v(:,2);
z = v(:,3);

plot3(x,y,z);
view(158, 14)
 

function uprime = lor2(t,u) 

global A B R 
uprime = zeros(3,1); 
uprime(1) = -A*u(1) + A*u(2); 
uprime(2) = R*u(1) - u(2) - u(1)*u(3); 
uprime(3) = -B*u(3) + u(1)*u(2);