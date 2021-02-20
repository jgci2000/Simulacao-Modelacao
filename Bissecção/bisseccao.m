clear all
close all
clc

%frep = 1;
w = 2 * pi;

g = 9.81;
y0 = 5;

t = 0:0.001:2;
res = zeros(2, length(t));
n = 1;

% for t = 0:0.001:2
%     res(1, n) = t;
%     res(2, n) = d(w, g, y0, t);
%     n = n + 1;
% end
% plot(res(1, :), res(2, :));
% axis([0 2 -1 11]);

t1 = 0;
t2 = 3;
delta = 1;
deltac = 0.0000001;
ta = t1;
tb = t2;
d1 = d(w, g, y0, ta);
d2 = d(w, g, y0, tb);

if d1 * d2 > 0
    disp("hello");
end

while delta > deltac
    
    t = (ta + tb) / 2;
    
    if d(w, g, y0, t) * d(w, g, y0, ta) < 0
        delta = abs(ta - t) / abs(t);
        tb = t;
    else
        delta = abs(tb - t) / abs(t);
        ta = t;
    end
    
end

disp(t);

for i = linspace(0, t, 100)
    
    if i == t
        
        for n = linspace(t, t + 1, 100)
            
            hold off
            
            vr = - w * sin(w * t);
            vb = g * t - vr;
            y0 = y0 - 0.5 * g * t * t;
            
            yr = cos(w * n);
            yb = y0 + vb * n - 0.5 * g * n * n;
            
            plot(0, yb, '.r', 'MarkerSize', 20);
            hold on
            plot([-1 1], [yr yr], '-k');
            axis([-2 2 -2 6]);
            pause(0.8)
        end
        
    end
    
    hold off
    yr = cos(w * i);
    yb = y0 - 0.5 * g * i * i;
    
    plot(0, yb, '.r', 'MarkerSize', 20);
    hold on
    plot([-1 1], [yr yr], '-k');
    axis([-2 2 -2 6]);
    pause(0.1)
    
end



function f = d(w, g, y0, t)
f = (y0 - 0.5 * g * t * t) - cos(w * t);
end