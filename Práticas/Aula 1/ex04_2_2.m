close all
clear all

r = 2;

hold on

for ang = 0:pi/50:2*pi
    
    plot(r * cos(ang) + 1, r * sin(ang), '--k')
    
    %xlim([-2.5 4.5])
    %ylim([-3 3])
    %pause(0.1)
end
hold off
