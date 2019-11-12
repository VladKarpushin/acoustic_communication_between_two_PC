% 2019-11-04
% equalizer
% Fs can be removed as input parameter
function [H_out] = equalizer(s_a, s_b, w, L_new)
s_a = s_a/max(s_a);
s_b = s_b/max(s_b);

x = 1:length(s_a);
x = x / length(x);

x_new = 1:L_new;
x_new = x_new / length(x_new);

H = (fft(s_a) ./ fft(s_b));
H = smooth(H, w);
H = smooth(H, w);
H_new = interp1(x, H, x_new, 'spline'); % H is complex
H_out = H_new;