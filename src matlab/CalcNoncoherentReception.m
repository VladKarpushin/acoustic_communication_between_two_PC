%2016-02-13 This function calculated Noncoherent reception
%реализуетс€ коррел€ционный приЄм по двум квадратурным составл€ющим.
function [SignAmpl SignalComplex] = CalcNoncoherentReception(SignR,T,F,Fs)
% input:
% 	SignR - received signal
%   T     - window of filter
%   Fs    - sample rate
%   F     - frequency of signal, 200<F<Fs/20000, [Hz]

% output:
% 	SignAmpl    - estimated amplitude.  вадратура огибающей
%   SignalComplex   - complex signal

Td = 2*pi/Fs;   %sampling interval
x=0:F*Td:length(SignR)*(F*Td);       %f0
x = x(1:length(SignR));
SignS = sin(x);          %sin signal
SignC = cos(x);          %cos signal

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
SignRC = SignR.*(SignC);

% x=1:length(SignR);
% figure, plot(x,SignR,'r',x,SignRS,'g',x,SignRC,'b');
% title('SignR (r) SignRS (g) and SignRC (b)');

SignRS_smooth = smooth(SignRS,T);   %Q(t) quadrature
SignRC_smooth = smooth(SignRC,T);   %I(t) in-phase

% x=1:length(SignR);
% figure, plot(x,SignR,'r',x,SignRS_smooth,'g',x,SignRC_smooth,'b');
% title('SignR (r) SignRS_smooth (g) and SignRC_smooth (b)');

%normalization (start)
%varSC   = var(SignS);
%varR    = var(SignR);
%varRSC  = sqrt(varR*varSC);
% SignRS_smooth = SignRS_smooth/varRSC; %in this case SignAmpl depends on SignR duration
% SignRC_smooth = SignRC_smooth/varRSC;

% SignRS_smooth = SignRS_smooth/std(SignS);   %in this case SignAmpl depends on SignR Amplitude
% SignRC_smooth = SignRC_smooth/std(SignS);
%normalization (stop)

SignAmpl = SignRS_smooth.^2+SignRC_smooth.^2;   %detected amplitude
%SignPhase = -atan(SignRS_smooth./SignRC_smooth);
SignalComplex = complex(SignRS_smooth, SignRC_smooth);

% x=1:length(SignR);
% figure, plot(x,SignR,'r',x,SignRS_smooth,'g',x,SignRC_smooth,'b',x,SignAmpl,'y',x,SignPhase,'k');
% title('SignR (r) SignRS_smooth (g) and SignRC_smooth (b) SignAmpl (y) SignPhase k');
% figure, plot(x,SignR,'r',x,SignRS_smooth,'g');
% title('SignR (r) SignRS_smooth (g)');
end
