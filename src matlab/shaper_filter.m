% 2016-12-13
% Function filters signal in order to reduce high frequencies 
% this is shaper function

function [signal_long_filtered] = shaper_filter(signal_long, samples, Fs)
% input:
% output:

%N=128;           %length of filter pulse response
N = samples * 4;
nfft = length(signal_long);
alpha = 8;      %For Gaussian window. The width of the window is inversely related to the value of ?. A larger value of ? produces a more narrow window
stdev = (N - 1) / (2 * alpha);    %For Gaussian window. Standard deviation, ?,

x = -pi:2 * pi/nfft:pi - pi / nfft;
windft_gausswinNew = exp(-1 / 2 * (x / (1 / stdev)) .^ 2)';
windft_gausswinNew = ifftshift(windft_gausswinNew);

N_rect = fix(nfft / samples);
windft_rectwin = zeros(nfft, 1);
windft_rectwin(1:N_rect) = 1;
windft_rectwin(nfft - N_rect:nfft) = 1;

windft = windft_gausswinNew; % windft can be replaced by G(f)
%windft = windft_rectwin;
signal_long_dft = fft(signal_long);
signal_long_filtered = real(ifft(signal_long_dft .* windft));

x=1:nfft;
x=x / Fs;
figure, plot(x, signal_long / max(signal_long), x, signal_long_filtered / max(signal_long_filtered));
xlabel('sec');
legend('signal_long', 'signal_long filterred');

signal_long_dft = signal_long_dft .* conj(signal_long_dft);
x = 1:nfft;
x = x / nfft * Fs;
figure, plot(x, signal_long_dft / max(signal_long_dft), x, abs(windft .* conj(windft)), x, windft_rectwin);
xlabel('Hz')
title('PSD of signal_long and windft and rectdft');
end