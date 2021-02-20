close all
clear all
hold on

r = 2;
c = [1; 0];
x = -1:0.1:3;
ang = 0:pi/50:2*pi
for x = -1:0.1:3
    
   %Gráfico 
   y1 = (sqrt((-(x-c(1)) .^ 2) + r .^2) + c(2));
   y2 = -(sqrt((-(x-c(1)) .^ 2) + r .^2) + c(2));
   plot(x, y1, '--k',  'lineWidth', 0.1)
   plot(x, y2, '--k', 'lineWidth', 0.1)
   xlim([-2.5 4.5])
   ylim([-3 3])
   
   %Ponto a mover-se
   
   py = sin(ang) * r;
   px = cos(ang) * r + 1;
   plot(px, py, 'ro')
   
   %yp1 = (sqrt((-(x-c(1)) .^ 2) + r .^2) + c(2));
   %yp2 = -(sqrt((-(x-c(1)) .^ 2) + r .^2) + c(2));
   %plot(x, yp1, 'ro')
   %plot(x, yp2, 'ro')
   %pause(0.01)
   
   
   
end

hold off