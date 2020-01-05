%2016-06-18 This function realises coherent reception algorithm (for BPSK
%only) by correlation integral calculation
%2016-07-06 added offset for PLL
%2016-11-05 added code for signal constellation

function [corr_integral signal_complex] = CalcCoherentReceptionNew3(z, T, F, Fs, PLL_offset_n)
% input:
% 	z - received signal (Rx signal)
%   T     - window of filter (quantity of samples per one symbol)
%   Fs    - sample rate
%   F     - frequency of signal, 200<F<Fs/20000, [Hz]
%   PLL_offset_n - %Phase Locked Loop(PLL) offset of sin wave. PLL_offset_n  = ceil(Samples/4) = Pi/4. PLL_offset_n = 0 means there is not offset
% output:
% 	output - cross correlation function (CCF) received signal z and sin wave. Another name is correlation integral


Td = 2 * pi / Fs;   %sampling interval
x = 0:F * Td:(length(z) + PLL_offset_n) * (F * Td);
x = x(1 + PLL_offset_n:length(z) + PLL_offset_n);
%SignS = sin(x);          %local oscullator, sin signal
%SignC = cos(x);          %local oscullator, cos signal

z_sin = z .* sin(x);        %coefficient of correlation with sin
z_cos = z .* cos(x);        %coefficient of correlation with cos

% z_sin = z .* (SignS);        %coefficient of correlation with sin
% z_cos = z .* (SignC);        %coefficient of correlation with cos
corr_integral = smooth(z_sin, T);      %CCF received signal and sin wave (correlation integral)
signal_complex = complex(corr_integral, smooth(z_cos, T));
%signal_complex = complex(z_sin, z_cos);
end
