%2016-01-12
%Sound transfer/receiver
%2016-03-10 zero errors during signal transmittion
%2016-04-04 transferred file from the same pc
%2016-04-26 finished
%2016-05-22 added duration of one symbol output/output
%2016-05-29 added time delay in a beginning of transmission (unit is bit)
%2016-06-05 added comments - Samples
%2016-06-18 auto-calculated threshold for OOK,
%2016-06-19 released BPSK with BER=0
%2016-07-03 started to develop "phase locking" (PLL - Phase Locked Loop)
%2016-07-06 realised PLL algorithm, tested for models
%2016-07-09 Successfully tested PLL with real signals 
%2016-07-22 Met problem with sync during >2KB transmitting

close all,clc,clear all;

%read a file start
%[signalInf_b errmsg] = file2signal('input\main13.m');
[signalInf_b errmsg] = file2signal('input\main13_2.0KB.m');
if length(errmsg) ~= 0
    disp('file2signal error');
    disp(errmsg);
    return;
end
%read a file stop

%******************************
%*******Transmitter************
%******************************
%Fs - Sampling rate in Hz. Valid values depend on the specific audio hardware installed. Typical values supported by most sound cards are 8000, 11025, 22050, 44100, 48000, and 96000 Hz.

%carrier signal forming(start)
Fs = 96000;     %sample rate

F = Fs/30;  %frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave
kt = 1;     %coefficient of duration of one symbol, kt/F = duration of one symbol
%nInfBits = 4*8*1024;   %number of information bits
nInfBits = 1000;   %number of information bits
%nInfBits = length(signalInf_b);

Td = 2*pi/Fs;   %sampling interval

SignBarker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
%SignBarker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13

delay = 0;       % time delay in a beginning of transmission (unit is bit)
nTotalBits = delay + 2*length(SignBarker) + nInfBits;
t = kt*nTotalBits/F;   %common transmit time

x = 0:F*Td:kt*nTotalBits*2*pi;
%uC0 = sin(x+9*pi/20);                    %carrier signal
uC0 = sin(x);                    %carrier signal
%uC0 = sin(x + pi/4);                    %carrier signal
%uC0 = -cos(x);                    %carrier signal
%figure,plot(uC0(1:30))

disp(['Sampling rate = ',num2str(Fs),' Hz']);
disp(['Freq of signal1 = ',num2str(F),' Hz']);
disp(['Duration of one symbol = ',num2str(1000*kt/F),' ms']);
disp(['BW = ',num2str(2/(kt/F)),' Hz']);
disp(['Number of total transmitted bits = ',num2str(nTotalBits),' [bits]']);
disp(['Number of information transmitted bits = ',num2str(nInfBits),' [bits]']);
disp(['Troughput = ',num2str(F/kt),' [bits/s]']);
disp(['Symbol rate = ',num2str(F/kt),' [Hz]']);
disp(['common transmit time = ',num2str(t),' [s]']);
%carrier signal forming(stop)

%modulation(start)
signalInf_b = 2*randi([0,1],nInfBits,1)-1; %information signal = noise

% signalInf_b = -1*ones(nInfBits,1);       %information signal = [1 -1 1 -1 1 -1]
% n=1:2:nInfBits;
% signalInf_b(n) = 1;
% 
% signalInf_b = ones(nInfBits,1);       %information signal = [-1 -1 -1 1 1 1]
% signalInf_b(1:fix(nInfBits/2)) = -1;

%adding sync marks (start)
signal = zeros(nTotalBits,1);           %sync signal + information signal + sync signal
signal(1+delay:length(SignBarker)+delay) = SignBarker;
signal(nTotalBits - length(SignBarker)+1:length(signal)) = -SignBarker;
signal(length(SignBarker)+1+delay:length(SignBarker)+delay + nInfBits) = signalInf_b;
%adding sync marks (stop)

Samples = kt*Fs/F;       %!!!! number of samples per one symbol
if abs(Samples - fix(Samples)) > 0                  %check Freq assignment error
    disp(['Error. Freq assignment error. Samples = kt*Fs/F = ',num2str(Samples)]);
    return;
end

u = uC0;
u = 0*u;        %u is transmitter signal
for n = 1+delay:nTotalBits
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples-1;
    if signal(n) == 1
        u(iA:iB) = uC0(iA:iB);   %BPSK
    else 
        u(iA:iB) = -uC0(iA:iB);   %BPSK
    end
end
%signal(1:5)
%modulation(stop)

x1 = 0:2*pi/100:kt*2*pi;
x2 = 0:F*Td:kt*2*pi;
figure,plot(x1,sin(x1),x2,sin(x2),'o');
%figure,plot(x1,uC0(x1),x2,uC0(x2),'o');
title('one symbol');

U_PSD = fft(u).*conj(fft(u));   %power spectrum density

x = 1:length(u);
x=x/Fs;
figure,plot(x,u);
xlabel('sec');
title('transmitted signal u');

x = 1:floor(length(U_PSD)/2);
figure,plot(x/t,U_PSD(x));
xlabel('Hz')
title('PSD of transmitted signal u');

nBits = 24;
sound(u,Fs,nBits);         %modulated signal

%return

%**************************************************
%***************Receiver***************************
%**************************************************

%Fs = 44100;% 44100; %96000 
% F = 5*Fs/100;  %frequency of signal, 200<F<Fs/2, [Hz]. even(1*Fs/100, 2*Fs/100, 4*Fs/100). F = 2(and 4)*Fs/100 - optimum
% kF = 4;         %f1 = kF*f0, kF=2(and 8) - optimum
% nInfBits = 1*8*1024;   %number of information bits
% SignBarker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
% %SignBarker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
tt = 1+kt*(2*length(SignBarker) + nInfBits)/F;   %common transmit time
nBits=24;
Samples = kt*Fs/F;       %!!!! number of samples per one symbol
if abs(Samples - fix(Samples)) > 0                  %check Freq assignment error
    disp(['Error. Freq assignment error. Samples = kt*Fs/F = ',num2str(Samples)]);
    return;
end

recObj = audiorecorder(Fs, nBits, 1);
get(recObj);
disp('Start recording.');
recordblocking(recObj, tt);
disp('End of Recording.');

% Store data in double-precision array.
z = getaudiodata(recObj)';      %received signal
%z = u;

SignBarkerLong = ones(length(SignBarker)*Samples,1);
for n = 1:length(SignBarker)
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    if SignBarker(n) == -1
        SignBarkerLong(iA:iB) = -1;
    end
end

% Plot the waveform.
x = 1:length(z);
x=x/Fs;
figure,plot(x,z);
xlabel('sec');
title('recorded signal z');

Z_PSD = fft(z).*conj(fft(z));   %power spectrum density
x = 1:floor(length(Z_PSD)/2);
figure,plot(x/tt,Z_PSD(x));
xlabel('Hz')
title('PSD of reveived signal z');

%*******PLL start ******
n = 2*ceil(Samples/4);          %number of PLL iterations n = Pi/4
PLL_offset_n = 0:n-1;           %PLL_offset_n - %sin offset (used for PLL). It is equal to PLL_offset_n  = ceil(Samples/4). PLL_offset_n = 0 means there is not offset

MaxAbsSignAmp = (-2)*ones(n,1);
BER =           (-2)*ones(n,1);
MaxSignSync =   (-2)*ones(n,1);
MinSignSync =   (-2)*ones(n,1);
delta = (-2)*ones(n,1);
threshold = 0;

for i = 1:n
    [SignAmp] = CalcCoherentReceptionNew(z,Samples,F,Fs,PLL_offset_n(i));   %coherent reception
    MaxAbsSignAmp(i) = max(abs(SignAmp));
    [EstSignal_b MaxSignSync(i) MinSignSync(i) Err delta(i)] = CalcSignalEstimationNew2(SignAmp,threshold, SignBarkerLong, Samples); %This function estimates information bits (information signal)
    if length(EstSignal_b) == length(signalInf_b)                  %check size
        BER(i) = mean(abs(EstSignal_b - signalInf_b)/2);   %The bit error rate (BER)
    else
        disp(['Error. Can not calculate BER, because of different array size. length(EstSignal_b) = ',num2str(length(EstSignal_b)), ',  length(signalInf_b) = ', num2str(length(signalInf_b))]);
    end

end
PLL_offset_vs_BER = [PLL_offset_n' MaxAbsSignAmp MaxSignSync MinSignSync MaxSignSync-MinSignSync BER delta abs(nInfBits*Samples - delta)];
m = -2;
i = -2;
[m i] = max(MaxSignSync-MinSignSync);
disp(['PLL_offset_n = ',num2str(PLL_offset_n(i))]);
disp(['MaxSignSync-MinSignSync = ',num2str(m)]);
disp(['BER = ',num2str(BER(i))]);
disp(['delta = ',num2str(delta(i))]);
[SignAmp] = CalcCoherentReceptionNew(z,Samples,F,Fs,PLL_offset_n(i));   %coherent reception
[EstSignal_b MaxSignSync MinSignSync Err] = CalcSignalEstimationNew2(SignAmp,threshold, SignBarkerLong, Samples); %This function estimates information bits (information signal)
x = 1:length(z);
%x=x/Fs;
figure,plot(x,SignAmp);
xlabel('sec');
title('SignAmp');

figure,plot(x,SignAmp,'r',x,z,'b');
xlabel('sec');
title('SignAmp (r) and z (b)');
%*******PLL stop ******

%****bit error calculation start*******
%The bit error rate (BER) is the number of bit errors per unit time. The bit error ratio (also BER) is the number of bit errors divided by the total number of transferred bits during a studied time interval. BER is a unitless performance measure, often expressed as a percentage.
% EstSignal_b
% signalInf_b
if length(EstSignal_b) ~= length(signalInf_b)                  %check Freq assignment error
    disp(['Error. Can not calculate BER, because of different array size. length(EstSignal_b) = ',num2str(length(EstSignal_b)), ',  length(signalInf_b) = ', num2str(length(signalInf_b))]);
    return;
end

BER = abs(EstSignal_b - signalInf_b)/2;   %The bit error rate (BER)
disp(['BER = ',num2str(mean(BER))]); 
disp(['Number of bit errors= ',num2str(sum(BER))]); 
%****bit error calculation stop*******

%write file start
[errmsg] = signal2file('output\output.txt', EstSignal_b);
%[errmsg] = signal2file(strcat('output\',char(datetime),'.txt'), EstSignal_b);
if length(errmsg) ~= 0
    disp('signal2file error');
    disp(errmsg);
    return;
end
%write file stop