% 2019-11-17
function [EstSignal_b, SignalContell, indexA, indexB, CorrIntegral, thr] = calc_ook_receiver_new(z,Samples,F,Fs, SignBarkerB1Long, SignBarkerB2Long, nInfBits, period, signalInf_b)

%****noncoherent reception start *******
[SignalComplex] = CalcNoncoherentReceptionNew(z,Samples,F,Fs);      %SignalComplex - complex signal
CorrIntegral = real(SignalComplex).^2+imag(SignalComplex).^2;       %detected amplitude (amplitude envelope quadrature)
%****noncoherent reception stop*******

%****information signal estimation (start)*******
threshold = 0:0.01:0.6;  %resolver threshold
n = length(threshold);
BER = (-2)*ones(n,1);
MaxSignSync = (-2)*ones(n,1);
MinSignSync = (-2)*ones(n,1);

for i = 1:n
    %[EstSignal_b MaxSignSync(i) MinSignSync(i) Err] = CalcSignalEstimationNew(CorrIntegral,threshold(i), SignBarkerLong, Samples); %This function estimates information bits (information signal)
    [EstSignal_b MaxSignSync(i) MinSignSync(i) Err] = CalcSignalEstimationNew4B1B2(CorrIntegral,threshold(i), SignBarkerB1Long,SignBarkerB2Long, Samples,nInfBits, period,SignalComplex); %This function estimates information bits (information signal)
    if length(EstSignal_b) == nInfBits                  %check Freq assignment error
        BER(i) = mean(abs(EstSignal_b - signalInf_b)/2);   %The bit error rate (BER)
    else
        disp(['Error. Can not calculate BER, because of different array size. length(EstSignal_b) = ',num2str(length(EstSignal_b)), ',  length(signalInf_b) = ', num2str(length(signalInf_b))]);
    end
end
%thr_vs_BER = [threshold' MaxSignSync MinSignSync BER];
%****information signal estimation (stop)*******

m = -2;
i = -2;
[m i] = max(MaxSignSync-MinSignSync);
disp(['threshold = ',num2str(threshold(i))]);
disp(['MaxSignSync-MinSignSync = ',num2str(m)]);
disp(['BER = ',num2str(BER(i))]);
[SignalComplex] = CalcNoncoherentReceptionNew(z,Samples,F,Fs);      %SignalComplex - complex signal
CorrIntegral = real(SignalComplex).^2+imag(SignalComplex).^2;       %detected amplitude (amplitude envelope quadrature)
[EstSignal_b a a a a a SignalContell indexA indexB] = CalcSignalEstimationNew4B1B2(CorrIntegral,threshold(i), SignBarkerB1Long,SignBarkerB2Long, Samples,nInfBits,period,SignalComplex); %This function estimates information bits (information signal)
thr = threshold(i);