% 2016-02-13 This function realise noncoherent reception
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

% figure;
% h = histogram(abs(z));
% title('z histogramm');
% hx = h.BinWidth*(1:h.NumBins);
% hy = h.Values;
% hn = round(h.NumBins*0.3);    %take 30%
% [m i] = max(hy(length(hy)-hn:length(hy)));
% EstMax = hx(length(hy)+i-hn-1);
% z = z / EstMax;
% figure,plot(z);
% title('z/EstMax');

%z = z/max(z);
z_sin = z .* sin(x);        %coefficient of correlation with sin
z_cos = z .* cos(x);        %coefficient of correlation with cos


% x=1:length(z);
% figure, plot(x,z,'r',x,z_sin,'g',x,z_cos,'b');
% title('z (r) z_sin (g) and z_cos (b)');

%zS_smooth = smooth(z_sin, T);   %Q(t) quadrature
%zC_smooth = smooth(z_cos, T);   %I(t) in-phase

% x=1:length(z);
% figure, plot(x,z,'r',x,zS_smooth,'g',x,zC_smooth,'b');
% title('z (r) zS_smooth (g) and zC_smooth (b)');

%normalization (start)
%varSC   = var(SignS);
%varR    = var(z);
%varRSC  = sqrt(varR*varSC);
% zS_smooth = zS_smooth/varRSC; %in this case SignAmpl depends on z duration
% zC_smooth = zC_smooth/varRSC;

% zS_smooth = zS_smooth/std(SignS);   %in this case SignAmpl depends on z Amplitude
% zC_smooth = zC_smooth/std(SignS);
%normalization (stop)

%SignAmpl = zS_smooth.^2+zC_smooth.^2;   %detected amplitude
%SignPhase = -atan(zS_smooth./zC_smooth);
%signal_complex = complex(zS_smooth, zC_smooth);
signal_complex = complex(smooth(z_cos, T), smooth(z_sin, T));

% x=1:length(z);
% figure, plot(x,z,'r',x,zS_smooth,'g',x,zC_smooth,'b',x,SignAmpl,'y',x,SignPhase,'k');
% title('z (r) zS_smooth (g) and zC_smooth (b) SignAmpl (y) SignPhase k');
% figure, plot(x,z,'r',x,zS_smooth,'g');
% title('z (r) zS_smooth (g)');
end
