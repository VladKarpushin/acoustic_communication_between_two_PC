% 2019-12-01
% equalizer input point

function [y_out] = equalizer_first(sign_x, sign_y, w, indexA)

sign_y_cut = sign_y(indexA-length(sign_x):indexA-1)';
H = equalizer(sign_x, sign_y_cut, w, length(sign_y));
y_out = real(ifft(fft(sign_y) .* H)); % should be (H)
y_out = y_out - mean(y_out);