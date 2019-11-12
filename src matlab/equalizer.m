% 2019-11-04
% equalizer
% Fs can be removed as input parameter
function [H_out] = equalizer(s_a, s_b, w, Fs, L_new)
s_a = s_a/max(s_a);
s_b = s_b/max(s_b);

x = 1:length(s_a);
x = x / length(x) * Fs;

x_new = 1:L_new;
x_new = x_new / length(x_new) * Fs;

H_a = abs(fft(s_a) ./ fft(s_b));
H_a = smooth(H_a, w);
H_a = smooth(H_a, w);
H_new = interp1(x, H_a, x_new, 'spline');
H_out = H_new;

% H_a = fft(s_a) ./ (fft(s_b) + 1);
% H_a = smooth(H_a, w);
% H_a = smooth(H_a, w);
% H_new = interp1(x, H_a, x_new, 'spline');
% H_out = H_new;