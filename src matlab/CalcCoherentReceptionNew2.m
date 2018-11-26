%2016-06-18 This function realises coherent reception algorithm (for BPSK
%only) by correlation integral calculation
%2016-07-06 added offset for PLL
function [output] = CalcCoherentReceptionNew2(SignR,T,F,Fs,PLL_offset_n)
% input:
% 	SignR - received signal (Rx signal)
%   T     - window of filter (quantity of samples per one symbol)
%   Fs    - sample rate
%   F     - frequency of signal, 200<F<Fs/20000, [Hz]
%   PLL_offset_n - %Phase Locked Loop(PLL) offset of sin wave. PLL_offset_n  = ceil(Samples/4) = Pi/4. PLL_offset_n = 0 means there is not offset
% output:
% 	output - cross correlation function (CCF) received signal SignR and sin wave. Another name is correlation integral


Td = 2*pi/Fs;   %sampling interval
x=0:F*Td:(length(SignR)+PLL_offset_n)*(F*Td);       %
x = x(1+PLL_offset_n:length(SignR)+PLL_offset_n);
SignS = sin(x);                 %local oscullator

SignRS = SignR.*(SignS);        %coefficient of correlation
output = smooth(SignRS,T);      %CCF received signal and sin wave (correlation integral)
end
