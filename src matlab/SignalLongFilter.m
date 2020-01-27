% 2016-12-13
% Function filters signal in order to reduce high frequencies 
% this is shaper function

function [SignalLong_filterred] = SignalLongFilter(SignalLong, Samples, Fs)
% input:
% output:

%N=128;           %length of filter pulse response
N = Samples*4;
nfft = length(SignalLong);
alpha = 8;      %For Gaussian window. The width of the window is inversely related to the value of ?. A larger value of ? produces a more narrow window
stdev = (N-1)/(2*alpha);    %For Gaussian window. Standard deviation, ?,

x = -pi:2*pi/nfft:pi-pi/nfft;
windft_gausswinNew = exp(-1/2*(x/(1/stdev)).^2)';
windft_gausswinNew = ifftshift(windft_gausswinNew);

N_rect = fix(nfft/Samples);
windft_rectwin = zeros(nfft,1);
windft_rectwin(1:N_rect) = 1;
windft_rectwin(nfft-N_rect:nfft) = 1;

windft = windft_gausswinNew;
%windft = windft_rectwin;
SignalLong_dft = fft(SignalLong);
SignalLong_filterred = real(ifft(SignalLong_dft.*windft));

x=1:nfft;
x=x/Fs;
figure, plot(x,SignalLong/max(SignalLong),x,SignalLong_filterred/max(SignalLong_filterred));
xlabel('sec');
legend('SignalLong','SignalLong filterred');

SignalLong_dft = SignalLong_dft.*conj(SignalLong_dft);
x = 1:nfft;
x = x/nfft*Fs;
figure,plot(x,SignalLong_dft/max(SignalLong_dft),x,abs(windft.*conj(windft)),x,windft_rectwin);
xlabel('Hz')
title('PSD of SignalLong and windft and rectdft');
end