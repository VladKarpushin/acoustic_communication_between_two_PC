%2016-05-12 OFDM transmitter. carrier-less OFDM over the entire 0-22KHz band
%2016-06-06 improved notations
%2016-06-07 realised OFDM receiver, added BER calculation
%2016-06-08 added white noise. Added jitter: algorithm is very sensitive to
%jitter
%2016-06-16 added OFDM signal BW = Sub-carrier BW*Nsc2_real, Number of real subcarriers Nsc2_real
%2016-10-17 realised inter-carrier aggregation (CA feature in LTE)
%2016-10-18 checked, we can use both: BPSK and OOK modulation
%2016-10-19 added resource grid
%2016-10-28 added signal constellation
%2016-11-01 added comments
%2016-11-07 added color for Signal Contellation
%2016-11-15 added SNR output [dB]
%2016-11-29 improved PSD figure
%future development - adding CP

close all,clc,clear all;


%******************************
%*******Input parameters*******
%******************************


Fs = 96000;                         %sample rate. Sampling rate in Hz. Valid values depend on the specific audio hardware installed. Typical values supported by most sound cards are 8000, 11025, 22050, 44100, 48000, and 96000 Hz.
%Fs = 48000;
Nsc = 128;                          %number of subcarriers
nInfBits = 1024*8*1;                %number of information bits, [bits]

%sc_list = [5,6,7,8];   %real subcarrier(frequency) list
%sc_list = [10];
sc_list = [5,6,7,8,9,11,12,13];   %real subcarrier(frequency) list

%SNR = 10e10;                        %signal to noise ratio, [times]
SNR = 1;                            %signal to noise ratio, [times]. How to convert it to dB? Very easy: SNR(dB) = 10*log10(P_signal/P_noise). For example, P_signal = P_noise, then SNR(dB) = 10*log10(1) = 0;
delay = 0;                          %signal delay, [samples]. It means OFDM symbol synchronization error.
threshold = 0;                      %resolver threshold of receiver (threshold  = 0 for FM, threshold = 0.5 for AM)

%******************************
%*******Transmitter************
%******************************
F = Fs/Nsc;                         %delta F, sub-carrier width
Nsc2 = Nsc/2;                       %half part of Nsc
Nsc2_real = length(sc_list);        %number of real subcarriers (frequencies)
Nsymb = nInfBits/Nsc2_real;         %number of symbols
if abs(Nsymb - fix(Nsymb)) > 0      %integer part verification
    disp(['Error. Nsymb calculation error. Nsymb = nInfBits/Nsc2_real = ',num2str(Nsymb)]);
    return;
end
t = Nsymb/F;                        %common transmission time, the same with length(signal_long)/Fs

signalInf_b = randi([0,1],nInfBits,1); %information signal = [rand]
%signalInf_b = zeros(nInfBits,1);       %information signal = [1 0 1 0 1 0]
%n=1:2:nInfBits;
%signalInf_b(n) = 1;

signal_s = zeros(Nsc2,Nsymb);

%serial to parallel convertation (S/P) (start)
i = 1;
for l = 1:Nsymb
    for j = 1:Nsc2_real
%        signal_s(sc_list(j),l) = signalInf_b(i);                       %AM OOK
        signal_s(sc_list(j),l) = 2*signalInf_b(i)-1;                    %FM BPSK sin/-sin
        i = i + 1;
    end
end
%serial to parallel convertation (S/P) (stop)

%adding symmetry part of spectrum (start)
%Properties of Fourier Transform: If a signal is even and real, then its spectrum is also real and even
%Properties of Fourier Transform: If a spectrum is even and real, then its signal is also real and even
signal_s = signal_s';
signal_s_flip = fliplr(signal_s);
signal_X = zeros(Nsc,Nsymb)';           %resource grid
signal_X(:,2:Nsc2+1) = signal_s;
signal_X(:,Nsc2+2:Nsc) = signal_s_flip(:,2:Nsc2);
%adding symmetry part of spectrum (stop)

signal = ifft(signal_X')';

%parallel to serial convertation (P/S) (start)
signal_long = zeros(Nsc*Nsymb,1);
ind = 1;
for i = 1:Nsymb
    signal_long(ind:ind+Nsc-1) = signal(i,:);
    ind = ind + Nsc;
end
signal_long = signal_long/std(signal_long); %signal with std=1
%parallel to serial convertation (P/S) (stop)

nBits = 24;
sound(signal_long,Fs,nBits);         %sound(y) sends audio signal y to the speaker 

%output results (start)
figure,imagesc(signal_X);        %to show resource grid
title('Resource grid');

x = 1:length(signal_long);
figure,plot(x/Fs,signal_long);
%figure,plot(signal_long);
xlabel('sec');
title('transmitted signal');

disp(['Sampling rate = ',num2str(Fs),' Hz']);
disp(['SNR = ',num2str(SNR),' [times]']);
disp(['SNR = ',num2str(10*log10(SNR)),' [dB]']);
disp(['Number of subcarriers Nsc = ',num2str(Nsc)]);
disp(['Number of real subcarriers Nsc2_real = ',num2str(Nsc2_real)]);
disp(['Number of symbols Nsymb = ',num2str(Nsymb)]);
disp(['Sub-carrier BW = ',num2str(F),' Hz']);
disp(['OFDM signal BW = Sub-carrier BW*Nsc2_real = ',num2str(F*Nsc2_real),' Hz']);
disp(['Duration of one symbol = ',num2str(1/F*1000),' ms']);
disp(['Common transmit time = ',num2str(t),' [s]']);
disp(['Troughput = ',num2str(nInfBits/t),' [bits/s]']);
disp(['Symbol rate = ',num2str(F),' [Hz]']);
%output results (stop)



%**************************************************
%***************Air channel************************
%**************************************************
%adding white noise and delay (start)
signal_noise = randn(size(signal_long))/SNR;  %SNR=3--->BER=0
signal_long = signal_long + signal_noise;
signal_long(delay+1:length(signal_long)) = signal_long(1:length(signal_long)-delay);
%adding white noise and delay (stop)

U_PSD = fft(signal_long).*conj(fft(signal_long));   %power spectrum density of signal
% x = 1:floor(length(U_PSD)/2);
% figure,plot((x-1)/t,U_PSD(x));
x = 1:length(U_PSD);
x = x/length(U_PSD)*Fs;
figure,plot(x,U_PSD);
xlabel('Hz')
title('PSD of signal+noise');

%**************************************************
%***************Receiver***************************
%**************************************************

EstSignal_b         = zeros(Nsc2_real*Nsymb,1);     %estimation of signalInf_b
SignalContell       = zeros(Nsc2_real*Nsymb,1);     %for signal constellation

signal_r = signal_long;                     %received signal

%serial to parallel convertation (S/P) (start)
signal_Y = zeros(Nsc,Nsymb);
signal_Y(1:length(signal_r)) = signal_r;
%serial to parallel convertation (S/P) (stop)

signal_Y = fft(signal_Y)';

%destroy symmetry part of signal
%Properties of Fourier Transform: If a signal is even and real, then its spectrum is also real and even
%Properties of Fourier Transform: If a spectrum is even and real, then its signal is also real and even

%parallel to serial convertation (P/S) (start)
i = 1;
for l = 1:Nsymb
    for j = 1:Nsc2_real
        SignalContell(i) = signal_Y(l,sc_list(j)+1)';
        i = i + 1;
    end
end
%parallel to serial convertation (P/S) (stop)

EstSignal_b = real(SignalContell);
EstSignal_b = EstSignal_b > max(EstSignal_b)*threshold;

%****bit error calculation start*******
%The bit error rate (BER) is the number of bit errors per unit time. The bit error ratio (also BER) is the number of bit errors divided by the total number of transferred bits during a studied time interval. BER is a unitless performance measure, often expressed as a percentage.
BER = abs(EstSignal_b - signalInf_b);   %The bit error rate (BER)
BERind = find(BER);                      %find - Find indices and values of nonzero elements
%****bit error calculation stop*******

%output results (start)
[signalInf_b(1:10) EstSignal_b(1:10)]
disp(['BER = ',num2str(mean(BER))]); 
disp(['Number of bit errors= ',num2str(sum(BER))]); 
figure, plot(real(signal_Y'));
title('Real part of received signal');
figure, plot(imag(signal_Y'));
title('Imag part of received signal');
% figure, plot(abs(signal_Y'),'k');
% title('Abs of received signal');

c = linspace(1,10,length(SignalContell));                   %from black to yellow
figure, scatter(real(SignalContell),imag(SignalContell),[],c);
ylim(xlim);
title('Signal Constellation');
xlabel('In phase');
ylabel('Quadrature');
%output results (stop)