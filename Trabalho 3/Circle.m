function Circle(x, y, r, c, lw)
% x and y are the coordinates of the center of the circle
% r is the radius of the circle
% 0.01 is the angle step, bigger values will draw the circle faster but
% you might notice imperfections (not very smooth)

% CIRCLE
ang = 0:0.01:2 * pi + 1;
xp = (r + lw) * cos(ang);
yp = (r + lw) * sin(ang);
plot(x + xp, y + yp, 'Color', c, 'LineWidth', lw);
end

