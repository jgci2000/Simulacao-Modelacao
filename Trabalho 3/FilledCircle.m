function FilledCircle(x, y, r, c)
% Filled Circle

rectangle('Position', [x - r, y - r, r * 2, r * 2], 'Curvature', [1 1], 'FaceColor', c);

end