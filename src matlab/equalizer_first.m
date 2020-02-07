% 2019-12-01
% equalizer
function [sign_y_out] = equalizer_first(sign_x, sign_y, w, index_a)

disp(['equalizer_first. index_a - length(sign_x) = ', num2str(index_a - length(sign_x))]);
disp(['equalizer_first. index_a - 1 = ', num2str(index_a - 1)]);
disp(['equalizer_first. len(sign_y) = ', num2str(length(sign_y))]);
sign_y_cut = sign_y(index_a - length(sign_x):index_a - 1)';
H = equalizer_second(sign_x, sign_y_cut, w, length(sign_y));
sign_y_out = real(ifft(fft(sign_y) .* H)); % should be (H)
sign_y_out = sign_y_out - mean(sign_y_out);
%figure, plot(abs(H));
%figure, plot(angle(H));
%figure, plot(sign_y_cut);
%index_a