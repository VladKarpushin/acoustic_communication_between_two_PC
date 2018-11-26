%2016-01-12
%Sound transfer/receiver
%2016-03-10 zero errors during signal transmittion
%2016-04-04 transferred file from the same pc
%2016-04-26 finished
%2016-05-22 added duration of one symbol output/output
%2016-05-29 added time delay in a beginning of transmission (unit is bit)
%2016-06-05 added comments - Samples
%2016-06-18 auto-calculated threshold
%2016-10-09 added new functions from BPSK version
%2016-10-11 realised sync algorithm. Now we can trasfer any number of bits
%2016-10-22 added netbitrate
%2016-11-20 updated modulation code
%2016-11-27 added code for signal constellation
%2016-11-30 improved PSD code
%2016-12-04 added code for SNR estimation
%2016-12-18 added filter. Successfully tested ultrasound
%2017-01-02 added spectrogram

tic
close all,clc,clear all;

%read a file start
filename = 'ones_1KB.m';
%[signalInf_b errmsg] = file2signal('input\ones_1KB.m');
[signalInf_b errmsg] = file2signal(strcat('input\',filename));
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
% Fs = 96000;     %sample rate
% F = Fs/35;  %frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave
Fs = 22050;     %sample rate
F = Fs/7;  %frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave
kt = 3;     %coefficient of duration of one symbol, kt/F = duration of one symbol
period = 1024*4;      %packet size
nInfBits = 1024*8*10;   %number of information bits
%nInfBits = length(signalInf_b);
Td = 2*pi/Fs;   %sampling interval
delay = 1000;       % time delay in a beginning of transmission (unit is bit)

%*****Barker codes set generation (start)*****
nSignBarkerB1 = 75;   %quantity of Barker codes in a set.
nSignBarkerB2 = 75;   %quantity of Barker codes in a set.
SignBarkerOneB1 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems. Barker codes have length at most 13 and have low correlation sidelobes
SignBarkerOneB2 = [1 1 1 -1 -1 -1 1 -1 -1 1 -1]'; %Barker code N=11. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems.
SignBarkerB1 = GetPeriodicBarkerCode(SignBarkerOneB1, nSignBarkerB1);
SignBarkerB2 = GetPeriodicBarkerCode(SignBarkerOneB2, nSignBarkerB2);
%*****Barker codes set generation (stop)*****

n = fix(nInfBits/period);
nTotalBits = delay + 2*length(SignBarkerB1) + n*length(SignBarkerB2) + nInfBits; %it needs to calculate
%nTotalBits = delay + 2*length(SignBarker) + nInfBits;
t = kt*nTotalBits/F;   %common transmit time

disp(['Sampling rate = ',num2str(Fs),' Hz']);
disp(['Freq of signal1 = ',num2str(F),' Hz']);
disp(['Duration of one symbol = ',num2str(1000*kt/F),' ms']);
disp(['BW = ',num2str(2/(kt/F)),' Hz']);
disp(['Number of total transmitted bits = ',num2str(nTotalBits),' [bits]']);
disp(['Number of information transmitted bits = ',num2str(nInfBits),' [bits]']);
disp(['Physical layer gross bitrate = ',num2str(F/kt),' [bits/s]']);
disp(['Net bitrate = ',num2str(fix(nInfBits/t)),' [bits/s]']);
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
%signal  = InsertSyncB1(signalInf_b,SignBarker, delay);
signal          = InsertSyncB2(signalInf_b,SignBarkerB2, period);
signal          = InsertSyncB1(signal,SignBarkerB1, delay);

%adding sync marks (stop)

Samples = kt*Fs/F;       %!!!! number of samples per one symbol
if abs(Samples - fix(Samples)) > 0                  %check Freq assignment error
    disp(['Error. Freq assignment error. Samples = kt*Fs/F = ',num2str(Samples)]);
    return;
end


%x = linspace(0,kt*nTotalBits*2*pi-(F*Td),nTotalBits*Samples);
x = 0:F*Td:(kt*nTotalBits*2*pi)-(F*Td);

SignalLong = (Short2Long(signal, Samples)+1)/2;     %OOK
%SignalLong = (Short2Long(signal, Samples));        %BPSK
SignalLong(1:delay*Samples) = 0;
SignalLong = SignalLongFilter(SignalLong, Samples, Fs);     %filtering
u = SignalLong.*sin(x)';

%signal(1:5)
%modulation(stop)

%air channel modeling (start)
u = u/std(u);
% SNR = 1000;
% signal_noise = randn(length(u),1)/SNR;
% u = u + signal_noise;
%air channel modeling (stop)


x1 = 0:2*pi/100:kt*2*pi;
x2 = 0:F*Td:kt*2*pi;
figure,plot(x1,sin(x1),x2,sin(x2),'o');
title('one symbol');

U_PSD = fft(u).*conj(fft(u));   %power spectrum density

x = 1:length(u);
x=x/Fs;
figure,plot(x,u);
xlabel('sec');
title('transmitted signal u');

% x = 1:floor(length(U_PSD)/2);
% figure,plot(x/t,U_PSD(x));
x = 1:length(U_PSD);
x = x/length(U_PSD)*Fs;
figure,plot(x,U_PSD);
xlabel('Hz')
title('PSD of transmitted signal u');

figure, spectrogram(u,400,100,[],Fs);
title('Transmitted signal spectrogram');

nBits = 24;
%sound([u.*0 u],Fs,nBits);         %modulated signal
sound(u,Fs,nBits);         %modulated signal

%return

%**************************************************
%***************Receiver***************************
%**************************************************

%info = audiodevinfo
%support = audiodevinfo(1,1,44100,16,1)
%Fs = 44100;% 44100; %96000 
% F = 5*Fs/100;  %frequency of signal, 200<F<Fs/2, [Hz]. even(1*Fs/100, 2*Fs/100, 4*Fs/100). F = 2(and 4)*Fs/100 - optimum
% kF = 4;         %f1 = kF*f0, kF=2(and 8) - optimum
% nInfBits = 1*8*1024;   %number of information bits
% SignBarker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
% %SignBarker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
%tt = 1+kt*(2*length(SignBarker) + nInfBits)/F;   %common transmit time
n = fix(nInfBits/period);
tt = 1+kt*(2*length(SignBarkerB1) + n*length(SignBarkerB2) + nInfBits)/F;   %common transmit time

%threshold = 0.12;   %resolver threshold
nBits=24;
Samples = kt*Fs/F;       %!!!! number of samples per one symbol
if abs(Samples - fix(Samples)) > 0                  %check Freq assignment error
    disp(['Error. Freq assignment error. Samples = kt*Fs/F = ',num2str(Samples)]);
    return;
end

recObj = audiorecorder(Fs, nBits, 1);
get(recObj);
% Record your voice for 5 seconds.
%recObj = audiorecorder;
disp('Start recording.');
recordblocking(recObj, tt+2);
disp('End of Recording.');

% Store data in double-precision array.
z = getaudiodata(recObj)';      %received signal
%z = u';

%SignBarkerLong = Short2Long(SignBarker, Samples);
SignBarkerB1Long    = Short2Long(SignBarkerB1, Samples);
SignBarkerB2Long    = Short2Long(SignBarkerB2, Samples);


% Plot the waveform.
x = 1:length(z);
x=x/Fs;
figure,plot(x,z);
xlabel('sec');
title('recorded signal z');

Z_PSD = fft(z).*conj(fft(z));   %power spectrum density
% x = 1:floor(length(Z_PSD)/2);
% figure,plot(x/tt,Z_PSD(x));
Z_PSD(length(Z_PSD)) = 0;
Z_PSD(1) = 0;
x = 1:length(Z_PSD);
x = x/length(Z_PSD)*Fs;
figure,plot(x,Z_PSD);

xlabel('Hz')
title('PSD of reveived signal z');

figure, spectrogram(z,400,100,[],Fs); % Compute the short-time Fourier transform. Divide the waveform into 400-sample segments with 100-sample overlap
title('Received signal spectrogram');


%****noncoherent reception start *******
[SignalComplex] = CalcNoncoherentReceptionNew(z,Samples,F,Fs);      %SignalComplex - complex signal
CorrIntegral = real(SignalComplex).^2+imag(SignalComplex).^2;       %detected amplitude (amplitude envelope quadrature)

% x = 1:length(z);
% x=x/Fs;
% figure,plot(x,CorrIntegral);
% xlabel('sec');
% title('SignAmp');
% 
% figure,plot(x,CorrIntegral,'r',x,z,'b');
% xlabel('sec');
% title('SignAmp (r) and z (b)');
%****noncoherent reception stop*******

%****information signal estimation (start)*******
%[EstSignal_b Err] = CalcSignalEstimation(SignAmp,threshold, SignBarkerLong, Samples, Fs); %This function estimates information bits (information signal)
threshold = 0:0.01:0.6;  %resolver threshold
n = length(threshold);
BER = (-2)*ones(n,1);
MaxSignSync = (-2)*ones(n,1);
MinSignSync = (-2)*ones(n,1);

for i = 1:n
    %[EstSignal_b MaxSignSync(i) MinSignSync(i) Err] = CalcSignalEstimationNew(CorrIntegral,threshold(i), SignBarkerLong, Samples); %This function estimates information bits (information signal)
    [EstSignal_b MaxSignSync(i) MinSignSync(i) Err] = CalcSignalEstimationNew4B1B2(CorrIntegral,threshold(i), SignBarkerB1Long,SignBarkerB2Long, Samples,nInfBits,period,SignalComplex); %This function estimates information bits (information signal)
    if length(EstSignal_b) == length(signalInf_b)                  %check Freq assignment error
        BER(i) = mean(abs(EstSignal_b - signalInf_b)/2);   %The bit error rate (BER)
    else
        disp(['Error. Can not calculate BER, because of different array size. length(EstSignal_b) = ',num2str(length(EstSignal_b)), ',  length(signalInf_b) = ', num2str(length(signalInf_b))]);
    end
end
thr_vs_BER = [threshold' MaxSignSync MinSignSync MaxSignSync-MinSignSync BER];
%****information signal estimation (stop)*******



%*******output result (start)*********
m = -2;
i = -2;
[m i] = max(MaxSignSync-MinSignSync);
disp(['threshold = ',num2str(threshold(i))]);
disp(['MaxSignSync-MinSignSync = ',num2str(m)]);
disp(['BER = ',num2str(BER(i))]);
[EstSignal_b a a a a a SignalContell indexA indexB] = CalcSignalEstimationNew4B1B2(CorrIntegral,threshold(i), SignBarkerB1Long,SignBarkerB2Long, Samples,nInfBits,period,SignalComplex); %This function estimates information bits (information signal)


indexA = indexA-length(SignBarkerB1Long);
indexB = indexB+length(SignBarkerB1Long);
if (indexA-4 > 1) && (indexB > 1) && (indexA < length(z)) && (indexB < length(z))  %SNR estimation 
    s = std(z(indexA:indexB));
    n = std(z(1:indexA-4));
    SNR_estimated = s/n;
    disp(['SNR estimated = ',num2str(round(SNR_estimated))]);
    disp(['SNR estimated = ',num2str(round(10*log10(SNR_estimated))), ' [dB]']);
end

x = 1:length(z);
x=x/Fs;
figure,plot(x,CorrIntegral);
xlabel('sec');
title('SignAmp');

figure,plot(x,CorrIntegral,'r',x,z,'b');
xlabel('sec');
title('SignAmp (r) and z (b)');

c = linspace(1,10,length(SignalContell));                   %from black to yellow
figure,scatter(real(SignalContell),imag(SignalContell),[],c);   %Create a scatter plot and vary the circle color.
hold on;
theta = linspace(0,2*pi);
r = sqrt(threshold(i));
x = r*cos(theta);
y = r*sin(theta);
plot(x,y);

%ylim(xlim);
axis equal; %Use the same length for the data units along each axis.
xlabel('In Phase');
ylabel('Quadrature');
title('Signal Constellation');
%*******output result (stop)*********


%****bit error calculation start*******
%The bit error rate (BER) is the number of bit errors per unit time. The bit error ratio (also BER) is the number of bit errors divided by the total number of transferred bits during a studied time interval. BER is a unitless performance measure, often expressed as a percentage.
% EstSignal_b
% signalInf_b
if length(EstSignal_b) ~= length(signalInf_b)                  %check Freq assignment error
    disp(['Error. Can not calculate BER, because of different array size. length(EstSignal_b) = ',num2str(length(EstSignal_b)), ',  length(signalInf_b) = ', num2str(length(signalInf_b))]);
    return;
end

% BER = abs(EstSignal_b - signalInf_b)/2;   %The bit error rate (BER)
% disp(['BER = ',num2str(mean(BER))]); 
% disp(['Number of bit errors= ',num2str(sum(BER))]); 
% %****bit error calculation stop*******

%write file start
%[errmsg] = signal2file('output\output.txt', EstSignal_b);
[errmsg] = signal2file(strcat('output\',filename), EstSignal_b);
%[errmsg] = signal2file(strcat('output\',char(datetime),'.txt'), EstSignal_b);
if length(errmsg) ~= 0
    disp('signal2file error');
    disp(errmsg);
    return;
end
%write file stop
toc