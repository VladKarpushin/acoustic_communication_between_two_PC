% 2019-11-04
% equalizer
function [H] = equalizer(s_a, s_b, Fs, w)

s_a = s_a/max(s_a);
s_b = s_b/max(s_b);

S_A_PSD = fft(s_a).*conj(fft(s_a));   %power spectrum density
S_B_PSD = fft(s_b).*conj(fft(s_b));   %power spectrum density

%w = round(w); % 75*14 = (Period of Barker code) * (number of smples per one symbol)
S_A_PSD = smooth(S_A_PSD, w);
S_B_PSD = smooth(S_B_PSD, w);

x = 1:length(s_a);
x = x / length(x) * Fs;

figure, plot(x, S_A_PSD);
xlabel('Hz')
title('PSD of sa');

figure, plot(x, S_B_PSD);
xlabel('Hz')
title('PSD of sb');
H = S_A_PSD ./ (S_B_PSD + 1);
figure, plot(x, H);
xlabel('Hz')
title('H');