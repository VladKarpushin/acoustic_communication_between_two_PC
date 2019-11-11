% 2019-11-04
% equalizer
function [H] = equalizer(s_a, s_b, w, Fs, L)

s_a = s_a/max(s_a);
s_b = s_b/max(s_b);

% S_A_PSD = fft(s_a).*conj(fft(s_a));   %power spectrum density
% S_B_PSD = fft(s_b).*conj(fft(s_b));   %power spectrum density

% S_A_PSD = smooth(S_A_PSD, w);
% S_B_PSD = smooth(S_B_PSD, w);

x = 1:length(s_a);
x = x / length(x) * Fs;

x_new = 1:L;
x_new = x_new / length(x_new) * Fs;

% figure, plot(x, S_A_PSD);
% xlabel('Hz')
% title('PSD of sa');
% 
% figure, plot(x, S_B_PSD);
% xlabel('Hz')
% title('PSD of sb');
% H = S_A_PSD ./ (S_B_PSD + 1);
% 
% figure, plot(x, H);
% xlabel('Hz')
% title('H');

H_a = abs(fft(s_a) ./ fft(s_b) + 1);
H_a = smooth(H_a, w)

figure, plot(x, smooth(H_a, w));
xlabel('Hz')
title('Ha');

H_new = interp1(x, H_a, x_new);
%H_new = interp1(x, H_a, x_new, 'spline');
figure, plot(x_new, H_new);
xlabel('Hz')
title('H new');

% figure, plot(x, smooth(H_ph, w));
% xlabel('Hz')
% title('Hph');

H = H_a;
