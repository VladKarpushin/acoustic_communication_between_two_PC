% 2019-02-02 This function calculates long complex signal

function [signal_complex] = calc_signal_complex(z, T, F, Fs, phi)
% input:
% 	z       - received signal (Rx signal)
%   T       - window of filter (quantity of samples per one BPSK/OOK symbol)
%   Fs      - sample rate
%   F       - carrier frequency of signal, 200<F<Fs/20000, [Hz]
%   phi     - phase of sin/cos
% output:
% 	signal_complex - coefficient of correlation with sin/cos
% 	wave. Another name is correlation integral or quadrature omponents of
% 	received signal

Td = 2 * pi / Fs;   % sampling interval
step = F * Td;
x = 0:step:(length(z) - 1) * step;

z = z / std(z);

z_sin = z .* sin(x + phi);        % coefficient of correlation with sin. Probably should be "-sin"
z_cos = z .* cos(x + phi);        % coefficient of correlation with cos

signal_complex = complex(smooth(z_cos, T), smooth(z_sin, T));
end
