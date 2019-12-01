% 2019-11-04
% equalizer
% Fs can be removed as input parameter
% Y = H * X, --> 1/H = X/Y

function [H_out] = equalizer_second(sign_x, sign_y, w, size_out)
sign_x = sign_x/max(sign_x);
sign_y = sign_y/max(sign_y);

H = fft(sign_x) ./ fft(sign_y);
H = smooth(H, w);
H = smooth(H, w);   % second time smoothing for better smoothing

x = 1:length(sign_x);
x = x / length(x);
x_out = 1:size_out;
x_out = x_out / length(x_out);
H_out = interp1(x, H, x_out, 'spline'); % H is complex