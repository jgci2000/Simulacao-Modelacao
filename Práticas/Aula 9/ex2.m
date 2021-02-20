clear all
close all
clc


mv = VideoReader('slow motion back flip.wmv');
nframes = int32(mv.Duration * mv.FrameRate);
frameno = 0;
n = 10;
t = 0;
c = 1;

for i = 1:(nframes / n)
    for j = 1:n
        mov = readFrame(mv);
        frameno = frameno + 1;
    end
    if frameno > 700 && frameno < 1100
        image(mov);
        title(strcat('Frame ', num2str(frameno), ' Ponto ', num2str(1)));
        [x1(i) y1(i)] = ginput(1);
        tempos(c) = t;
        t = t + 0.1;
        c = c + 1;
    end
end

plot(tempos, y1(71: end), '-r');
a = polyfit(tempos, y1(71:end), 2);
y = a(3) + a(2) * tempos + a(1) * tempos .^ 2;
hold on
plot(tempos, y, '- b')


