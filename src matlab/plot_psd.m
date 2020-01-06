% 2019-12-08
% Power spectrum density (PSD) display function

function plot_psd(z, Fs, str_x_label, str_title)

z_spectrum = fft(z);
%Z_PSD = fft(z) .* conj(fft(z));   % power spectrum density
Z_PSD = z_spectrum .* conj(z_spectrum);   % power spectrum density
Z_PSD(length(Z_PSD)) = 0;
Z_PSD(1) = 0;
x = 0:length(Z_PSD) - 1;
x = x / length(Z_PSD) * Fs;
figure, plot(x, Z_PSD);
xlabel(str_x_label)
title(str_title);
