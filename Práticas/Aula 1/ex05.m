close all
clear all
axis equal


for x = 0:0.1:5
    
    quiver(0 , 0, 5, 5, 0, '.k')
    hold on
    plot(x, x, 'ro')
    pause(0.1)
    hold off
    
end

for x = 5:-0.1:0
    
    quiver(5 , 5, -5, -5, 0, '.k')
    hold on
    plot(x, x, 'ro')
    pause(0.1)
    hold off
    
end

