% 2016-06-18 This function realises coherent reception algorithm (for BPSK
% only) by correlation integral calculation
% 2016-07-06 added offset for PLL
% 2016-11-05 added code for signal constellation
% 2020-01-08 added phi
% can be replaced by calc_signal_complex()
function [signal_complex] = calc_coherent_reception(z, T, F, Fs, phi)
% input:
% 	z - received signal (Rx signal)
%   T     - window of filter (quantity of samples per one symbol)
%   Fs    - sample rate
%   F     - frequency of signal, 200<F<Fs/20000, [Hz]
%   phi     - phase of sin
% output:
% 	output - cross correlation function (CCF) received signal z and sin wave. Another name is correlation integral

Td = 2 * pi / Fs;   %sampling interval
step = F * Td;
x = 0:step:(length(z) - 1) * step;

z = z / std(z);

z_sin = z .* sin(x + phi);        %coefficient of correlation with sin
z_cos = z .* cos(x + phi);        %coefficient of correlation with cos

signal_complex = complex(smooth(z_cos, T), smooth(z_sin, T));
end
