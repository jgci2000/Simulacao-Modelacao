clear all
close all
clc

A1 = [0, 0; 0, 0.16];
b1 = [0; 0];

A2 = [0.2 -0.26; 0.23 0.22];
b2 = [0; 1.6];

A3 = [-0.15 0.28; 0.26 0.24];
b3 = [0; 0.44];

A4 = [0.85 0.04; -0.04 0.85];
b4 = [0; 1.6];

v0 = [0; 0];

A = [A1, A2, A3, A4];
b = [b1, b2, b3, b4];

res = zeros(2, 1000);
res(:, 1) = [0; 0];

for i = 2:1000
    
    a_idx = randi(8);
    if mod(a_idx, 2) == 0
        a_idx = a_idx - 1;
    end
    b_idx = randi(4);
    res(:, i) = (A(:, a_idx:a_idx + 1)) * res(:, i - 1) + b(:, b_idx);
    
end

plot(res(1, :), res(2, :));

