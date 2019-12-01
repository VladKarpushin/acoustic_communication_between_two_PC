% 2019-12-01
% equalizer input point

function [sign_y_out] = equalizer_first(sign_x, sign_y, w, indexA)

sign_y_cut = sign_y(indexA-length(sign_x):indexA-1)';
H = equalizer_second(sign_x, sign_y_cut, w, length(sign_y));
sign_y_out = real(ifft(fft(sign_y) .* H)); % should be (H)
sign_y_out = sign_y_out - mean(sign_y_out);