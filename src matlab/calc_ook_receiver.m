% 2019-11-17
function [thr_vs_BER] = calc_ook_receiver(z,Samples,F,Fs, SignBarkerB1Long, SignBarkerB2Long, nInfBits, period)

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
    if length(EstSignal_b) == length(signalInf_b)                  %check Freq assignment error
        BER(i) = mean(abs(EstSignal_b - signalInf_b)/2);   %The bit error rate (BER)
    else
        disp(['Error. Can not calculate BER, because of different array size. length(EstSignal_b) = ',num2str(length(EstSignal_b)), ',  length(signalInf_b) = ', num2str(length(signalInf_b))]);
    end
end
thr_vs_BER = [threshold' MaxSignSync MinSignSync MaxSignSync-MinSignSync BER];
%****information signal estimation (stop)*******