%2016-01-12
%Sound transfer/receiver
%2016-03-10 zero errors during signal transmittion
%2016-04-04 transferred file from the same pc
%2016-04-26 finished
%2016-05-22 added duration of one symbol output/output
%2016-05-29 added time delay in a beginning of transmission (unit is bit)
%2016-06-05 added comments - samples
%2016-06-18 auto-calculated threshold for OOK,
%2016-06-19 released BPSK with BER=0
%2016-07-03 started to develop "phase locking" (PLL - Phase Locked Loop)
%2016-07-06 realised PLL algorithm, tested for models
%2016-07-09 Successfully tested PLL with real signals 
%2016-07-22 Met problem with sync during >2KB transmitting
%2016-07-24 std(VKF) output, improved PLL criteria
%2016-07-27 added comments
%2016-08-13 added GetPeriodicBarkerCode, InsertSyncB2 and InsertSyncB1
%2016-08-21 added SignBarkerLong = Short2Long(sign_barker, samples);
%2016-08-25 added Long2Short
%2016-11-05 added code for signal constellation
%2016-11-06 added color for phase
%2016-11-20 updated modulation code
%2016-11-30 improved PSD code
%2016-12-13 added tx filter code
%2016-12-20 added code for SNR estimation
%2017-01-02 added spectrogram
%2019-11-04 start to equalizer development

close all, clc, clear all;

%read a file start
%[signal_inf_bits errmsg] = file2signal('input\main13.m');
%[signal_inf_bits, errmsg] = file2signal('..\input\main13_2KB.m');
filename = 'ones_1KB.m';
[signal_inf_bits errmsg] = file2signal(strcat('..\input\',filename));

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
%Fs = 96000;
Fs = 22050;     %sample rate
F = Fs/7;  %frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave
%F = Fs/5;  %frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave
kt = 2;     %coefficient of duration of one symbol, kt/F = duration of one symbol
n_inf_bits = 1024 * 4 * 1;      % number of information bits
%n_inf_bits = length(signal_inf_bits);
Td = 2 * pi / Fs;   % sampling interval
delay = 1000;       % time delay in a beginning of transmission (unit is bit)
SNR = 10;        %signal to noise ratio

%*****Barker codes set generation (start)*****
n_sign_barker = 75;   %quantity of Barker codes in a set.
sign_barker_one = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems. Barker codes have length at most 13 and have low correlation sidelobes
sign_barker = GetPeriodicBarkerCode(sign_barker_one, n_sign_barker);
%*****Barker codes set generation (stop)*****

n_total_bits = delay + 2*length(sign_barker) + n_inf_bits;
show_sign_para(kt, F, Fs, n_total_bits, n_inf_bits);

%modulation(start)
signal_inf_bits = 2 * randi([0, 1], n_inf_bits, 1) - 1; % model of information signal is noise

% signal_inf_bits = -1*ones(n_inf_bits,1);       %information signal = [1 -1 1 -1 1 -1]
% n=1:2:n_inf_bits;
% signal_inf_bits(n) = 1;
% 
% signal_inf_bits = ones(n_inf_bits,1);       %information signal = [-1 -1 -1 1 1 1]
% signal_inf_bits(1:fix(n_inf_bits/2)) = -1;

%adding sync marks (start)
signal  = InsertSyncB1(signal_inf_bits,sign_barker, delay);
%adding sync marks (stop)

samples = kt * Fs / F;       %!!!! number of samples per one symbol
if abs(samples - fix(samples)) > 0                  % check Freq assignment error
    disp(['Error. Freq assignment error. samples = kt*Fs/F = ', num2str(samples)]);
    return;
end

x = 0:F * Td:(kt * n_total_bits * 2 * pi) - (F * Td);
%x = linspace(0,kt*n_total_bits*2*pi-(F*Td),n_total_bits*samples);
%signal_long = (Short2Long(signal, samples)+1)/2;     %OOK
signal_long = (Short2Long(signal, samples));        %BPSK
signal_long(1:delay * samples) = 0;
%signal_long = SignalLongFilter(signal_long, samples, Fs);     %filtering

u = signal_long.*sin(x)';

%modulation(stop)


%air channel modeling (start)
u = u/std(u);
% signal_noise = randn(length(u),1)/SNR;
% u = u + signal_noise;
%air channel modeling (stop)


x1 = 0:2 * pi / 100:kt * 2 * pi;
x2 = 0:F * Td:kt * 2 * pi;
figure, plot(x1, sin(x1), x2, sin(x2), 'o');
title('one symbol');

plot_time(u, Fs, 'sec', 'transmitted signal u')
plot_psd(u, Fs, 'Hz', 'PSD of transmitted signal u');

figure, spectrogram(u, 400, 100, [], Fs);
title('Transmitted signal spectrogram');

nBits = 24;
sound(u, Fs, nBits);         %modulated signal

%**************************************************
%***************Receiver***************************
%**************************************************

%Fs = 44100;% 44100; %96000 
% F = 5*Fs/100;  %frequency of signal, 200<F<Fs/2, [Hz]. even(1*Fs/100, 2*Fs/100, 4*Fs/100). F = 2(and 4)*Fs/100 - optimum
% kF = 4;         %f1 = kF*f0, kF=2(and 8) - optimum
% n_inf_bits = 1*8*1024;   %number of information bits
tt = 1+kt*(2*length(sign_barker) + n_inf_bits)/F;   %common transmit time
nBits=24;
samples = kt*Fs/F;       %!!!! number of samples per one symbol
if abs(samples - fix(samples)) > 0                  %check Freq assignment error
    disp(['Error. Freq assignment error. samples = kt*Fs/F = ',num2str(samples)]);
    return;
end

recObj = audiorecorder(Fs, nBits, 1);
get(recObj);
disp('Start recording.');
recordblocking(recObj, tt);
disp('End of Recording.');

% Store data in double-precision array.
z = getaudiodata(recObj)';      %received signal
%z = u';

% SignBarkerLong = ones(length(sign_barker)*samples,1);
% for n = 1:length(sign_barker)
%     iA = 1 + (n-1)*samples;
%     iB = iA+samples-1;
%     if sign_barker(n) == -1
%         SignBarkerLong(iA:iB) = -1;
%     end
% end

SignBarkerLong = Short2Long(sign_barker, samples);

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


%*******PLL start ******
n = 2*ceil(samples/4);          %number of PLL iterations n = Pi/4
PLL_offset_n = 0:n-1;           %PLL_offset_n - %sin offset (used for PLL). It is equal to PLL_offset_n  = ceil(samples/4). PLL_offset_n = 0 means there is not offset

MaxAbsCorrIntegral   = (-2)*ones(n,1);   %max(CCF received signal and sin wave), max(correlation integral)
BER             = (-2)*ones(n,1);   %BER is bit error rate
MaxSignSync     = (-2)*ones(n,1);   %MaxSignSync is max(CCF)
MinSignSync     = (-2)*ones(n,1);   %MinSignSync is min(CCF)
delta           = (-2)*ones(n,1);   %delta is difference between index(MinSignSync) and index(MaxSignSync)
StdSignSync     = (-2)*ones(n,1);   %StdSignSync is std(CCF)
threshold = 0;                      %resolver threshold. Should be zero for BPSK

for i = 1:n
    [CorrIntegral tmp] = CalcCoherentReceptionNew3(z,samples,F,Fs,PLL_offset_n(i));   %coherent reception
    MaxAbsCorrIntegral(i) = max(abs(CorrIntegral));
    [EstSignal_b MaxSignSync(i) MinSignSync(i) Err delta(i) StdSignSync(i)] = CalcSignalEstimationNew4(CorrIntegral,threshold, SignBarkerLong, samples, tmp); %This function estimates information bits (information signal)
    if length(EstSignal_b) == length(signal_inf_bits)                  %check size
        BER(i) = mean(abs(EstSignal_b - signal_inf_bits)/2);   %The bit error rate (BER) calculation
    else
        disp(['Error. Can not calculate BER, because of different array size. length(EstSignal_b) = ',num2str(length(EstSignal_b)), ',  length(signal_inf_bits) = ', num2str(length(signal_inf_bits))]);
    end

end
ErrSyst = n_inf_bits*samples - delta; %systematic error between n_inf_bits*samples and delta
PLL_offset_vs_BER = [PLL_offset_n' MaxAbsCorrIntegral MaxSignSync MinSignSync MaxSignSync-MinSignSync BER delta ErrSyst StdSignSync];
%*******PLL stop ******

%*******output result (start)*********
m = -2;
i = -2;
[m i] = max(MaxSignSync-MinSignSync);   %PLL criterion1: max(CCF) - min(CCF) = max
if abs(ErrSyst(i))>samples/2      
    [m i] = min(abs(ErrSyst));          %PLL criterion2: systematic error = min
end

disp(['PLL_offset_n = ',num2str(PLL_offset_n(i))]);
disp(['BER = ',num2str(BER(i))]);
disp(['MaxSignSync - MinSignSync = ',num2str(MaxSignSync(i) - MinSignSync(i))]);
disp(['delta = ',num2str(delta(i))]);
disp(['ErrSyst (systematic error) = ',num2str(ErrSyst(i))]);
disp(['StdSignSync = ',num2str(StdSignSync(i))]);

[CorrIntegral SignalComplex] = CalcCoherentReceptionNew3(z,samples,F,Fs,PLL_offset_n(i));   %coherent reception
[EstSignal_b MaxSignSync MinSignSync Err delta StdSignSync SignalContell indexA indexB] = CalcSignalEstimationNew4(CorrIntegral,threshold, SignBarkerLong, samples, SignalComplex); %This function estimates information bits (information signal)

% equalizer start()
%sign_x = SignalLongFilter(SignBarkerLong, samples, Fs);     %filtering
sign_x = SignBarkerLong;
x = 0:F*Td:(kt*n_total_bits*2*pi)-(F*Td);
sign_x = sign_x.*sin(x(1:length(sign_x)))';

z_new = equalizer_first(sign_x, z, 3 * n_sign_barker, indexA);

Z_new_PSD = fft(z_new).*conj(fft(z_new));   %power spectrum density
Z_new_PSD(1) = 0;
x = 1:length(z_new);
x = x/length(z_new)*Fs;
figure, plot(x, Z_new_PSD);
xlabel('Hz')
title('PSD of equalized z');

[CorrIntegral SignalComplex] = CalcCoherentReceptionNew3(z_new,samples,F,Fs,PLL_offset_n(i));   %coherent reception
[EstSignal_b MaxSignSync MinSignSync Err delta StdSignSync SignalContell indexA indexB] = CalcSignalEstimationNew4(CorrIntegral,threshold, SignBarkerLong, samples, SignalComplex); %This function estimates information bits (information signal)
BER_eq = mean(abs(EstSignal_b - signal_inf_bits)/2);   %The bit error rate (BER) calculation
disp(['BER_eq = ', num2str(BER_eq)]);
% equalizer stop()


indexA = indexA-length(SignBarkerLong);
indexB = indexB+length(SignBarkerLong);
if (indexA-4 > 1) && (indexB > 1) && (indexA < length(z)) && (indexB < length(z))  %SNR estimation 
    s = std(z(indexA:indexB));
    n = std(z(1:indexA-4));
    SNR_estimated = s/n;
    disp(['SNR estimated = ',num2str(round(SNR_estimated))]);
    disp(['SNR estimated = ',num2str(round(10*log10(SNR_estimated))), ' [dB]']);
end

x = 1:length(z);
%x=x/Fs;
figure,plot(x,CorrIntegral);
xlabel('sec');
title('CorrIntegral');

figure,plot(x,CorrIntegral,'r',x,z,'b');
xlabel('sec');
title('CorrIntegral (r) and z (b)');

c = linspace(1,10,length(SignalContell));                   %from black to yellow
figure,scatter(real(SignalContell),imag(SignalContell),[],c);   %Create a scatter plot and vary the circle color.
%figure,scatter(real(SignalContell),imag(SignalContell),25,c,'filled');
ylim(xlim);
xlabel('In Phase');
ylabel('Quadrature');
title('Signal Constellation');
%*******output result (stop)*********

%****bit error calculation start*******
%The bit error rate (BER) is the number of bit errors per unit time. The bit error ratio (also BER) is the number of bit errors divided by the total number of transferred bits during a studied time interval. BER is a unitless performance measure, often expressed as a percentage.
% EstSignal_b
% signal_inf_bits
% if length(EstSignal_b) ~= length(signal_inf_bits)                  %check Freq assignment error
%     disp(['Error. Can not calculate BER, because of different array size. length(EstSignal_b) = ',num2str(length(EstSignal_b)), ',  length(signal_inf_bits) = ', num2str(length(signal_inf_bits))]);
%     return;
% end
% 
% BER = abs(EstSignal_b - signal_inf_bits)/2;   %The bit error rate (BER)
% disp(['BER = ',num2str(mean(BER))]); 
% disp(['Number of bit errors= ',num2str(sum(BER))]); 
% %****bit error calculation stop*******
% 

%write file (start)
[errmsg] = signal2file('output\output.txt', EstSignal_b);
%[errmsg] = signal2file(strcat('output\',char(datetime),'.txt'), EstSignal_b);
if length(errmsg) ~= 0
    disp('signal2file error');
    disp(errmsg);
    return;
end
%write file (stop)