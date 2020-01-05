%function [est_signal_b Err] = CalcSignalEstimation(CorrIntegrall,threshold)
%This function estimates information bits (information signal)
%without plots
%for BPSK only
%2016-11-05 added code for signal constellation
%2016-12-20 added iAB1 and iBB1 export for SNR calculation
function [est_signal_b  MaxSignSync MinSignSync Err delta StdSignSync signal_constel iA iB] = CalcSignalEstimationNew4(corr_integral, threshold, sign_barker_long, samples, signal_complex)
% input:
% 	corr_integral    - cross correlation function (CCF) received signal SignR and sin wave. Another name is correlation integral
% 	threshold       - resolver threshold. Should be zero for BPSK
%   sign_barker_long  -    sync signal(long)
%   samples         - quantity of samples per one symbol
%   signal_complex   - complex signal
% output:
% 	est_signal_b     - estimated information bits
%   Err             - error information
%   MaxSignSync     - MaxSignSync is max(CCF)
%   MinSignSync     - MinSignSync is min(CCF)
%   delta           - delta is difference between index(MinSignSync) and index(MaxSignSync)
%   StdSignSync     - StdSignSync is std(CCF) between two mainlobes
%   signal_constel   - signal constellation
%   iAB1            - index of first symbol. For SNR estimation
%   iBB1            - index of last symbol. For SNR estimation

Err = 0;
MaxSignSync = 0;
MinSignSync = 0;
est_signal_b = 0;
signal_constel = 0;
iA = 0;
iB = 0;

EstSignal = zeros(length(corr_integral),1);
EstSignal = (2 * (corr_integral > threshold))-1;    %resolver
% x = 1:length(EstSignal);
% x=x/Fs;
% figure,plot(x,EstSignal);
% title('EstSignal');

%****syncronization start*******
[SignSync Err] = CalcCCF_FFT(EstSignal,sign_barker_long,0);
% figure,plot(SignSync);
% title('SignSync');

[MaxSignSync,ImaxSignSync] = max(SignSync);  %largest element index
[MinSignSync,IminSignSync] = min(SignSync);  %smalles element index
% disp('Sync signal information');
% disp(['MaxSignSync = ',num2str(MaxSignSync),', ImaxSignSync = ',num2str(ImaxSignSync)]);
% disp(['MinSignSync = ',num2str(MinSignSync),', IminSignSync = ',num2str(IminSignSync)]);

iA = ImaxSignSync + length(sign_barker_long);
iB = IminSignSync-1;
if iA > iB                  
    iA = IminSignSync + length(sign_barker_long);
    iB = ImaxSignSync-1;
    EstSignal = -EstSignal;
    signal_complex = -signal_complex;
%     disp('Error. Sync error. ImaxSignSync > IminSignSync ');
%     Err = 1;
%     return;
end
StdSignSync = std(SignSync(iA:iB-length(sign_barker_long)));  %std(CCF) between two mainlobes
% SignSyncAdd = zeros(length(SignSync),1);
% SignSyncAdd(iA) = 1;
% SignSyncAdd(iB-length(sign_barker_long)) = 1;
% x = 1:length(SignSync);
%figure, plotyy(x,SignSync,x,SignSyncAdd)

est_signal_b = Long2Short(EstSignal(iA:iB), samples);
signal_constel = Long2Short(signal_complex(iA:iB), samples);
delta = iB - iA;

% nMax = round((iB-iA)/samples);
% iA = iA + fix(samples/2);
% est_signal_b =  zeros(nMax,1);
% i = 1;
% for n = iA:samples:iB
%     est_signal_b(i) = EstSignal(n);
%     i = i + 1;
% end
%****syncronization stop*******

if abs(length(est_signal_b)/8 - fix(length(est_signal_b)/8)) > 0                  %checking if length of est_signal_b is multiple with 8
    disp(['Error. abs(length(est_signal_b)/8 - fix(length(est_signal_b)/8)) > 0. length(est_signal_b) = ',num2str(length(est_signal_b))]);
    Err = 1;
    return;
end
% x=1:length(EstSignalDecimation);
% figure,stem(EstSignalDecimation);
% title('EstSignalDecimation');
% figure,plot(x,EstSignalDecimation,x,EstSignal);
% title('EstSignalDecimation + EstSignal');
end
