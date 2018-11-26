%2017-01-29
%PLL freq and phse adjustment investigation

close all,clc,clear all;

Fs = 22050;     %sample rate
F = Fs/100;       %frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave
Size = Fs;      %size of signal
dF = 10;        %[Hz] - frequency offset between Local oscillator and received signal

Td = 2*pi/(Fs+dF);  %sampling interval of received signal
x = 0:F*Td:Size*F*Td;
Sign_Rx = sin(x);    %received signal         


Td = 2*pi/Fs;   %sampling interval in Local oscillator
x_local = 0:F*Td:Size*F*Td;
Sign_local = sin(x_local);          

figure,plot(x_local,Sign_Rx,'r',x_local,Sign_local,'b');
legend('SignRx','SignLocal');