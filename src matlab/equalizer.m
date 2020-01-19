% 2019-12-01
% equalizer input point

function [z_out] = equalizer(z, ind_a, sign_barker_b1_long, n_sign_barker, samples, Fs, F)

sign_x = SignalLongFilter(sign_barker_b1_long, samples, Fs);     %filtering
%sign_x = sign_barker_long;

Td = 2 * pi / Fs;   % sampling interval

x = (0:length(sign_x) - 1) * F * Td;
sign_x = sign_x .* cos(x)';
%x = 0:F * Td:(kt * n_total_bits * 2 * pi) - (F * Td);
%sign_x = sign_x .* cos(x(1:length(sign_x)))';

z_out = equalizer_first(sign_x, z, 3 * n_sign_barker, ind_a);
%plot_psd(z_out, Fs, 'Hz', 'PSD of equalized received z');
