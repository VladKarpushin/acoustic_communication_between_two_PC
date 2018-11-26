%2016-06-18 This function realises coherent reception algorithm (for BPSK only)
%2016-07-06 added offset for PLL
function [SignRS_smooth] = CalcCoherentReceptionNew(SignR,T,F,Fs,PLL_offset_n)
% input:
% 	SignR - received signal (Rx signal)
%   T     - window of filter (quantity of samples per one symbol)
%   Fs    - sample rate
%   F     - frequency of signal, 200<F<Fs/20000, [Hz]
%   PLL_offset_n - %Phase Locked Loop(PLL) offset of sin wave. PLL_offset_n  = ceil(Samples/4) = Pi/4. PLL_offset_n = 0 means there is not offset
% output:
% 	SignRS_smooth    - cross correlation function (CCF) received signal
% 	SignR and sin  wave. Another name of this function is Q(t) quadrature

Td = 2*pi/Fs;   %sampling interval
x=0:F*Td:(length(SignR)+PLL_offset_n)*(F*Td);       %
x = x(1+PLL_offset_n:length(SignR)+PLL_offset_n);
SignS = sin(x);          %sin signal

%we search peak in 30% signal (from ceil) for solving PARP problem
% [hyy,edges] = histcounts(abs(SignR)); %partitions the X values into bins, and returns the count in each bin, as well as the bin edges. The histcounts function uses an automatic binning algorithm that returns bins with a uniform width, chosen to cover the range of elements in X and reveal the underlying shape of the distribution.
% hxx = edges(2:length(edges));   %maybe this string can be removed
% hnn = round(length(hxx)*0.3);    %take 30%
% [m i] = max(hyy(length(hyy)-hnn:length(hyy)));
% EstMax = hxx(length(hyy)+i-hnn-1);

% h = histogram(abs(SignR));
% hx = h.BinWidth*(1:h.NumBins);
% hy = h.Values;
% hn = round(h.NumBins*0.3);    %take 30%
% [m i] = max(hy(length(hy)-hn:length(hy)));
% EstMax = hx(length(hy)+i-hn-1);

%SignR = SignR/EstMax;
% figure,plot(SignR);
% title('SignR/EstMax');

%SignR = SignR/max(SignR);
SignRS = SignR.*(SignS);
SignRS_smooth = smooth(SignRS,T);   %CCF
end
