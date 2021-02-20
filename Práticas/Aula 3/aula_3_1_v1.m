clear all
close all

A = 2;
lambda_2 = 1:0.25:6;
fi = 0;

for lambda = 1:0.25:6
    for x = 0:12
        x_2 = x / lambda;
        y = (A * sin(2 * pi * x_2 + fi));
        %plot(f(x, lambda, fi, A), x, '-b')
        plot(y, '-b')
        xlim = [0, 12];
        hold on
    end
end

%function y = f(x, lambda, fi, A)
%y = (A * sin(2*pi*(x/lambda) + fi));
%end


