%2016-06-18 This function calculated Coherent reception
%реализуется когерентный приёмник
function [SignRS_smooth] = CalcCoherentReception(SignR,T,F,Fs,PLL_n)
% input:
% 	SignR - received signal
%   T     - window of filter
%   Fs    - sample rate
%   F     - frequency of signal, 200<F<Fs/20000, [Hz]
%   PLL_n - %sin offset (used for PLL). It is equal to PLL_n  =
%   ceil(Samples/4). PLL_n = 0 means there is not offset

% output:
% 	SignRS_smooth    - 

Td = 2*pi/Fs;   %sampling interval
x=0:F*Td:(length(SignR))*(F*Td);       %f0
%x=0:F*Td:(length(SignR)+PLL_n)*(F*Td);       %f0
x = x(1:length(SignR));
%x = x(1+PLL_n:length(SignR)+PLL_n);
SignS = sin(x);          %sin signal
%SignC = cos(x);          %cos signal

figure;
h = histogram(abs(SignR));
title('z histogramm');
hx = h.BinWidth*(1:h.NumBins);
hy = h.Values;
hn = round(h.NumBins*0.3);    %take 30%
[m i] = max(hy(length(hy)-hn:length(hy)));
EstMax = hx(length(hy)+i-hn-1);
SignR = SignR/EstMax;
figure,plot(SignR);
title('SignR/EstMax');

%SignR = SignR/max(SignR);
SignRS = SignR.*(SignS);
%SignRC = SignR.*(SignC);

% x=1:length(SignR);
% figure, plot(x,SignR,'r',x,SignRS,'g',x,SignRC,'b');
% title('SignR (r) SignRS (g) and SignRC (b)');

SignRS_smooth = smooth(SignRS,T);   %Q(t) quadrature
% x=1:length(SignR);
% figure, plot(x,SignR,'r',x,SignRS_smooth,'g');
% title('SignR (r) SignRS_smooth (g)');
end
