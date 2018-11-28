%2016-05-12 OFDM transmitter

close all,clc,clear all;

nInfBits = 8;
%signalInf_b = randi([0,1],nInfBits,1); %information signal = noise
signalInf_b = [0 1 2 3 10 3 2 1]; %information signal = noise
signalInf_b = [0 1 1 1 10 1 1 1]; %information signal = noise

U_PSD = fft(signalInf_b).*conj(fft(signalInf_b));   %power spectrum density
ifft(signalInf_b)
fft(signalInf_b)

figure,plot(U_PSD);
% xlabel('sec');
% title('transmitted signal u');
