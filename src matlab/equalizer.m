% 2019-11-04
% equalizer
function [H] = equalizer(s_a, s_b, Fs)

s_a = s_a/max(s_a);
s_b = s_b/max(s_b);

S_A_PSD = fft(s_a).*conj(fft(s_a));   %power spectrum density
S_B_PSD = fft(s_b).*conj(fft(s_b));   %power spectrum density

x = 1:length(s_a);
x = x / length(x) * Fs;

figure, plot(x, S_A_PSD);
xlabel('Hz')
title('PSD of sa');

figure, plot(x, S_B_PSD);
xlabel('Hz')
title('PSD of sb');
H = S_A_PSD ./ S_B_PSD;
figure, plot(x, H);
xlabel('Hz')
title('H');