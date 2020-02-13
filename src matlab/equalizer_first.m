% 2019-12-01
% equalizer
function [sign_y_out] = equalizer_first(sign_x, sign_y, w, index_a)

sign_y_out = sign_y;
ind_a = index_a - length(sign_x);
ind_b = index_a - 1;

if length(sign_y) < ind_b - ind_a
    disp(['equalizer_first. ind_a = ', num2str(ind_a)]);
    disp(['equalizer_first. ind_b = ', num2str(ind_b)]);
    disp(['equalizer_first. len(sign_y) = ', num2str(length(sign_y))]);
    disp(['length(sign_y) < inb_b - ind_a = ', num2str(ind_b - ind_a)]);
    return;
end

sign_y_cut = sign_y(ind_a:ind_b)';
H = equalizer_second(sign_x, sign_y_cut, w, length(sign_y));
sign_y_out = real(ifft(fft(sign_y) .* H)); % should be (H)
sign_y_out = sign_y_out - mean(sign_y_out);
%figure, plot(abs(H));
%figure, plot(angle(H));
%figure, plot(sign_y_cut);
%index_a