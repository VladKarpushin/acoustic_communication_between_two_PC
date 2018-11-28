%2016-05-12 OFDM transmitter
%2016-06-06 improved notations
%2016-06-07 realised OFDM receiver, added BER calculation
%2016-06-08 added white noise. Added jitter: algorithm is very sensitive to
%jitter
%2016-06-16 added OFDM signal BW = Sub-carrier BW*Nsc2_real, Number of real subcarriers Nsc2_real
%2016-10-15 sucessfully transferred from work PC to home PC

close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

%Fs = 96000;             %sample rate. Sampling rate in Hz. Valid values depend on the specific audio hardware installed. Typical values supported by most sound cards are 8000, 11025, 22050, 44100, 48000, and 96000 Hz.
Fs = 48000;
Nsc = 128;              %number of subcarriers
nInfBits = 1024;        %number of information bits
F = Fs/Nsc;             %delta F ,sub-carrier width
Nsc2 = Nsc/2;
% i_sc_start = 13;         %1 - min
% i_sc_stop  = Nsc2-20;      %Nsc2 - max
i_sc_start = 13;         %1 - min
i_sc_stop  = i_sc_start+3;      %Nsc2 - max
% i_sc_start = 1;         %1 - min
% i_sc_stop  = Nsc2;      %Nsc2 - max
Nsc2_real = i_sc_stop - i_sc_start + 1; %number of real subcarriers (Nsc2_real<=Nsc2)
    
%Nsymb = nInfBits/Nsc2;   %number of symbols
Nsymb = nInfBits/Nsc2_real;   %number of symbols
if abs(Nsymb - fix(Nsymb)) > 0                  %check Freq assignment error
    disp(['Error. Nsymb calculation error. Nsymb = nInfBits/Nsc2_real = ',num2str(Nsymb)]);
    return;
end


t = Nsymb/F;   %common transmit time, the same with length(signal_long)/Fs
%t = length(signal_long)/Fs;   %common transmit time, the same with Nsymb/F

%signal_X = randi([0,1],nInfBits,1); %information signal = noise
%signal_X = [0 1 2 3 10 3 2 1]; %information signal = noise

%signal_s = randi([0,1],nInfBits,1);
%signal_s = [1 2 3 4 5 6 7 8];
%signal_s = 1:nInfBits;

signalInf_b = randi([0,1],nInfBits,1);
%signalInf_b = 2*randi([0,1],nInfBits,1)-1;
% signalInf_b = zeros(nInfBits,1);
% signalInf_b(3) = 1;

signal_s = zeros(Nsc2,Nsymb);  
%signal_s(1:nInfBits) = signalInf_b;
i = 1;
for l = 1:Nsymb
    for k = i_sc_start:i_sc_stop
        %signal_s(k,l) = signalInf_b(i);        %AM OOK
        signal_s(k,l) = 2*signalInf_b(i)-1;     %BPSK sin/-sin
        i = i + 1;
    end
end

%signal_s(1:nInfBits) = 1:nInfBits;
% signal_s(10,1) = 1;     %resource element
% signal_s(20,2) = 1;
% signal_s(30,3) = 1;
% signal_s(40,4) = 1;

%making symmetry part of spectrum (start)
%Properties of Fourier Transform: If a signal is even and real, then its spectrum is also real and even
%Properties of Fourier Transform: If a spectrum is even and real, then its signal is also real and even
signal_s = signal_s';
signal_s_flip = fliplr(signal_s);
signal_X = zeros(Nsc,Nsymb)';   %resource block/grid
signal_X(:,2:Nsc2+1) = signal_s;
signal_X(:,Nsc2+2:Nsc) = signal_s_flip(:,2:Nsc2);
%making symmetry part of spectrum (start)

% signal_X(1+k) = 1;
% signal_X(Nsc+1-k) = 1;

signal = ifft(signal_X')';
% signal(:,1:2) = 0;
% signal(:,Nsc) = 0;
%signal(:,1) = 0;    %to destroy DC
signal_long = zeros(Nsc*Nsymb,1);
ind = 1;
for i = 1:Nsymb
    signal_long(ind:ind+Nsc-1) = signal(i,:);
    ind = ind + Nsc;
end
signal_long = signal_long/std(signal_long); %to normalize std

%adding white noise and jitter (start)
% signal_noise = randn(size(signal_long))/3;  %SNR=3--->BER=0
% signal_long = signal_long + signal_noise;
% jitter = 2;
% signal_long(jitter:length(signal_long)) = signal_long(1:length(signal_long)-jitter+1);
%adding white noise and jitter(stop)


%x = 1:length(signal_long);
%x=x/Fs;
%figure,plot(x,signal_long);
figure,plot(signal_long);
%xlabel('sec');
title('transmitted signal u');

U_PSD = fft(signal_long).*conj(fft(signal_long));   %power spectrum density
x = 1:floor(length(U_PSD)/2);
figure,plot((x-1)/t,U_PSD(x));
xlabel('Hz')
title('PSD of transmitted signal u');

disp(['Sampling rate = ',num2str(Fs),' Hz']);
disp(['Number of subcarriers Nsc = ',num2str(Nsc)]);
disp(['Number of real subcarriers Nsc2_real = ',num2str(Nsc2_real)]);
disp(['Number of symbols Nsymb = ',num2str(Nsymb)]);
disp(['Sub-carrier BW = ',num2str(F),' Hz']);
disp(['OFDM signal BW = Sub-carrier BW*Nsc2_real = ',num2str(F*Nsc2_real),' Hz']);
disp(['Duration of one symbol = ',num2str(1/F*1000),' ms']);
disp(['Common transmit time = ',num2str(t),' [s]']);
disp(['Troughput = ',num2str(nInfBits/t),' [bits/s]']);
disp(['Symbol rate = ',num2str(F),' [Hz]']);

%**************************************************
%***************Receiver***************************
%**************************************************

%threshold = 0.5;
threshold = 0.5;
signal_r = signal_long;     %received signal
signal_Y = zeros(Nsc,Nsymb);
signal_Y(1:length(signal_r)) = signal_r;
signal_Y = fft(signal_Y)';
%signal_Y(:,1) = 0;      %to destroy DC
figure, plot(real(signal_Y'));
title('Real part of received signal');
figure, plot(imag(signal_Y'));
title('Imag part of received signal');

signal_Y_real = real(signal_Y')';
EstSignal_b = zeros(Nsc2_real*Nsymb,1);  %estimation of signalInf_b
%destroy symmetry part of signal
%Properties of Fourier Transform: If a signal is even and real, then its spectrum is also real and even
%Properties of Fourier Transform: If a spectrum is even and real, then its signal is also real and even
ind = 1;
for i = 1:Nsymb
%     EstSignal_b(ind:ind+Nsc2-1) = signal_Y_real(i,2:Nsc2+1);
%     ind = ind + Nsc2;
    EstSignal_b(ind:ind+Nsc2_real-1) = signal_Y_real(i,i_sc_start+1:i_sc_stop+1);
    ind = ind + Nsc2_real;
end
EstSignal_b = EstSignal_b > max(EstSignal_b)*threshold;
%EstSignal_b = EstSignal_b > max(EstSignal_b)/10;   %the best for noise
%EstSignal_b = EstSignal_b > 0;

%****bit error calculation start*******
%The bit error rate (BER) is the number of bit errors per unit time. The bit error ratio (also BER) is the number of bit errors divided by the total number of transferred bits during a studied time interval. BER is a unitless performance measure, often expressed as a percentage.
[signalInf_b(1:10) EstSignal_b(1:10)]
BER = abs(EstSignal_b - signalInf_b);   %The bit error rate (BER)
BERind = find(BER)                      %find - Find indices and values of nonzero elements
disp(['BER = ',num2str(mean(BER))]); 
disp(['Number of bit errors= ',num2str(sum(BER))]); 
%****bit error calculation stop*******
