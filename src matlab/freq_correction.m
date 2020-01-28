 % 2020-01-03
% estimation of freq offset befween tx and rx, Fs estimation, frequency
% correction

function [est_F, est_Fs] = freq_correction(z, ind_a, synch_burst_size, freq_burst_size, Fs, F)

est_F = F;
est_Fs = Fs;
if freq_burst_size <= 10
    disp(['freq_burst_size <= 10']);
    return;
end

disp(['freq_correction debug. a = ', num2str(1 + ind_a - synch_burst_size - freq_burst_size)]);
disp(['freq_correction debug. b = ', num2str(ind_a - synch_burst_size)]);
z = z(1 + ind_a - synch_burst_size - freq_burst_size:ind_a - synch_burst_size)';

n = length(z) * 1000;
%n = 2^nextpow2(length(z) * 1000);
z_spectrum = fft(z, n);
%Z_PSD = fft(z, n) .* conj(fft(z, n));   % power spectrum density
Z_PSD = z_spectrum .* conj(z_spectrum);   % power spectrum density

%Z_PSD = fft(z) .* conj(fft(z));   % power spectrum density
Z_PSD(length(Z_PSD)) = 0; % can be replaced to "Z_PSD(end) = 0;"
%Z_PSD(1) = 0;
Z_PSD(1:fix(0.01 * length(Z_PSD))) = 0; % in order to reduce zero power
x = 0:length(Z_PSD) - 1;
x = x / length(Z_PSD) * Fs;
x = x(1:length(Z_PSD) / 2);
Z_PSD = Z_PSD(1:length(Z_PSD) / 2);

[m i] = max(Z_PSD);

est_phi = angle(z_spectrum(i));
est_F = x(i);
delta = 1 - est_F / F;
est_Fs = Fs * (1 + delta);

disp(['estimated phi = ', num2str(est_phi  * 180 / pi), ' [degrees]']);
disp(['estimated F = ', num2str(est_F), ' Hz']);
disp(['estimated F precision = ', num2str(x(2)), ' Hz']);
disp(['est_F - F = ', num2str(est_F - F), ' Hz']);
disp(['est_Fs - Fs = ', num2str(est_Fs - Fs), ' Hz']);

figure, plot(x, Z_PSD);
xlabel('Hz')
title('PLL signal PSD');
