%function [EstSignal_b Err] = CalcSignalEstimation(SignAmpl,threshold)
%This function estimated information bits (information signal)
%without plots
%for OOK
function [EstSignal_b MaxSignSync MinSignSync Err] = CalcSignalEstimationNew(SignAmp, threshold, SignBarkerLong, Samples)
% input:
% 	SignAmp - amplitude envelope quadrature
% 	threshold - resolver threshold
%   SignBarkerLong -     
%   Samples
%   Fs - 
% output:
% 	EstSignal_b - estimated information bits
%   Err - error information

Err = 0;
MaxSignSync = 0;
MinSignSync = 0;
EstSignal_b = 0;

EstSignal = zeros(length(SignAmp),1);
EstSignal = (2*(SignAmp > threshold))-1;    %for passive pause
% x = 1:length(EstSignal);
% x=x/Fs;
% figure,plot(x,EstSignal);
% title('EstSignal');

%****syncronization start*******
%[SignSync Err] = VKPCalcVKP_FFT(EstSignal,SignBarkerLong,0);
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
if iA > iB                  %check sync error
    disp('Error. Sync error. ImaxSignSync > IminSignSync ');
    Err = 1;
    return;
end

if MaxSignSync < 0                  %check sync error
    disp('Error. Sync error. MaxSignSync < 0');
    Err = 1;
    return;
end


EstSignal_b     = Long2Short(EstSignal(iA:iB),Samples);
% 
% nMax = round((iB-iA)/Samples);
% iA = iA + fix(Samples/2);
% EstSignalDecimation = zeros(length(EstSignal),1);
% EstSignal_b =  zeros(nMax,1);
% i = 1;
% for n = iA:Samples:iB
%     EstSignalDecimation(n) = EstSignal(n);
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
