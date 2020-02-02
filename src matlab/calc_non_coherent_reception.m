% 2016-02-13 This function realises noncoherent reception
% реализуетс€ коррел€ционный приЄм по двум квадратурным составл€ющим.




function [signal_complex] = calc_non_coherent_reception(z, T, F, Fs)
% input:
% 	z - received signal
%   T     - window of filter
%   Fs    - sample rate
%   F     - frequency of signal, 200<F<Fs/20000, [Hz]

% output:
%   signal_complex   - complex signal

Td = 2 * pi / Fs;   %sampling interval
step = F * Td;
x = 0:step:(length(z) - 1) * step;

z = z / std(z);

z_sin = z .* sin(x);        %coefficient of correlation with sin
z_cos = z .* cos(x);        %coefficient of correlation with cos

signal_complex = complex(smooth(z_cos, T), smooth(z_sin, T));
end
