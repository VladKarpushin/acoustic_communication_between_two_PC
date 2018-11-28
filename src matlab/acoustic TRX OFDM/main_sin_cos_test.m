%this script show spectr of sin and cos
%2016-10-26

close all,clc,clear all;
N=128;
s = zeros(N,1);
s(2) = 1;
s(N) = 1;
sx = ifft(s);
sX = fft(sx);

figure, plot(sx/max(sx));


x=1:128;
y=cos(2*pi*(x-1)/N);
figure,plot(y);

Y = fft(y);

PSD = fft(y).*conj(fft(y));   %power spectrum density

figure, plot(PSD);

[sx/max(sx) y']

figure, plot(x,sx/max(sx),x,y);