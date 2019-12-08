% 2019-12-08
% simple plot function

function plot_time(z, Fs, str_x_label, str_title)

x = 1:length(z);
x = x / Fs;
figure, plot(x, z);
xlabel(str_x_label)
title(str_title);
