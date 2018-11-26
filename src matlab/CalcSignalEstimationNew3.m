%function [EstSignal_b Err] = CalcSignalEstimation(CorrIntegrall,threshold)
%This function estimates information bits (information signal)
%without plots
%for BPSK only
%2016-11-05 added code for signal constellation
function [EstSignal_b  MaxSignSync MinSignSync Err delta StdSignSync SignalContell] = CalcSignalEstimationNew3(CorrIntegral, threshold, SignBarkerLong, Samples, SignalComplex)
% input:
% 	CorrIntegral    - cross correlation function (CCF) received signal SignR and sin wave. Another name is correlation integral
% 	threshold       - resolver threshold. Should be zero for BPSK
%   SignBarkerLong  -    sync signal(long)
%   Samples         - quantity of samples per one symbol
%   SignalComplex   - complex signal
% output:
% 	EstSignal_b     - estimated information bits
%   Err             - error information
%   MaxSignSync     - MaxSignSync is max(CCF)
%   MinSignSync     - MinSignSync is min(CCF)
%   delta           - delta is difference between index(MinSignSync) and index(MaxSignSync)
%   StdSignSync     - StdSignSync is std(CCF) between two mainlobes
%   SignalContell   - signal constellation

Err = 0;
MaxSignSync = 0;
MinSignSync = 0;
EstSignal_b = 0;
SignalContell = 0;

EstSignal = zeros(length(CorrIntegral),1);
EstSignal = (2*(CorrIntegral > threshold))-1;    %resolver
% x = 1:length(EstSignal);
% x=x/Fs;
% figure,plot(x,EstSignal);
% title('EstSignal');

%****syncronization start*******
[SignSync Err] = CalcCCF_FFT(EstSignal,SignBarkerLong,0);
% figure,plot(SignSync);
% title('SignSync');

[MaxSignSync,ImaxSignSync] = max(SignSync);  %largest element index
[MinSignSync,IminSignSync] = min(SignSync);  %smalles element index
% disp('Sync signal information');
% disp(['MaxSignSync = ',num2str(MaxSignSync),', ImaxSignSync = ',num2str(ImaxSignSync)]);
% disp(['MinSignSync = ',num2str(MinSignSync),', IminSignSync = ',num2str(IminSignSync)]);

iA = ImaxSignSync + length(SignBarkerLong);
iB = IminSignSync-1;
if iA > iB                  
    iA = IminSignSync + length(SignBarkerLong);
    iB = ImaxSignSync-1;
    EstSignal = -EstSignal;
    SignalComplex = -SignalComplex;
%     disp('Error. Sync error. ImaxSignSync > IminSignSync ');
%     Err = 1;
%     return;
end
StdSignSync = std(SignSync(iA:iB-length(SignBarkerLong)));  %std(CCF) between two mainlobes
% SignSyncAdd = zeros(length(SignSync),1);
% SignSyncAdd(iA) = 1;
% SignSyncAdd(iB-length(SignBarkerLong)) = 1;
% x = 1:length(SignSync);
%figure, plotyy(x,SignSync,x,SignSyncAdd)

EstSignal_b = Long2Short(EstSignal(iA:iB),Samples);
SignalContell = Long2Short(SignalComplex(iA:iB),Samples);
delta = iB-iA;

% nMax = round((iB-iA)/Samples);
% iA = iA + fix(Samples/2);
% EstSignal_b =  zeros(nMax,1);
% i = 1;
% for n = iA:Samples:iB
%     EstSignal_b(i) = EstSignal(n);
%     i = i + 1;
% end
%****syncronization stop*******

if abs(length(EstSignal_b)/8 - fix(length(EstSignal_b)/8)) > 0                  %checking if length of EstSignal_b is multiple with 8
    disp(['Error. abs(length(EstSignal_b)/8 - fix(length(EstSignal_b)/8)) > 0. length(EstSignal_b) = ',num2str(length(EstSignal_b))]);
    Err = 1;
    return;
end
% x=1:length(EstSignalDecimation);
% figure,stem(EstSignalDecimation);
% title('EstSignalDecimation');
% figure,plot(x,EstSignalDecimation,x,EstSignal);
% title('EstSignalDecimation + EstSignal');
end
