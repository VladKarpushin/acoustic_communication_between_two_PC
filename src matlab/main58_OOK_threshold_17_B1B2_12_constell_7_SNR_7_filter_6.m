%2016-01-12
%Sound transfer/receiver
%2016-03-10 zero errors during signal transmittion
%2016-04-04 transferred file from the same pc
%2016-04-26 finished
%2016-05-22 added duration of one symbol output/output
%2016-05-29 added time delay in a beginning of transmission (unit is bit)
%2016-06-05 added comments - samples
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
%2019-11-16 start to equalizer development

tic
close all, clc, clear all;

%read a file start
file_name = 'ones_1KB.m';
%[signal_inf_bits errmsg] = file2signal('input\ones_1KB.m');
[signal_inf_bits errmsg] = file2signal(strcat('..\input\', file_name));
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
% F = Fs/5;  %frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave

Fs = 22050;         % sample rate
F = Fs / 7;         % frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave
kt = 2;             % coefficient of duration of one symbol, kt/F = duration of one symbol
period = 1024 * 4 * 1;          % packet size
n_inf_bits = 1024 * 4 * 2;      % number of information bits

%n_inf_bits = length(signal_inf_bits);
Td = 2 * pi / Fs;   % sampling interval
delay = 1000;       % time delay in a beginning of transmission (unit is bit)

%*****Barker codes set generation (start)*****
n_sign_barker_b1 = 75;   % quantity of Barker codes in a set.
n_sign_barker_b2 = 75;   % quantity of Barker codes in a set.
sign_barker_one_b1 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; % Barker code N=13. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems. Barker codes have length at most 13 and have low correlation sidelobes
sign_barker_one_b2 = [1 1 1 -1 -1 -1 1 -1 -1 1 -1]'; % Barker code N=11. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems.
sign_barker_b1 = GetPeriodicBarkerCode(sign_barker_one_b1, n_sign_barker_b1);
sign_barker_b2 = GetPeriodicBarkerCode(sign_barker_one_b2, n_sign_barker_b2);
%*****Barker codes set generation (stop)*****

n = fix(n_inf_bits / period);
n_total_bits = delay + 2 * length(sign_barker_b1) + n * length(sign_barker_b2) + n_inf_bits; % total bits plus delay
%n_total_bits = delay + 2*length(SignBarker) + n_inf_bits;

show_sign_para(kt, F, Fs, n_total_bits, n_inf_bits);
%carrier signal forming(stop)

% modulation(start)
signal_inf_bits = 2 * randi([0, 1], n_inf_bits, 1) - 1; % model of information signal is noise

% signal_inf_bits = -1*ones(n_inf_bits,1);       % information signal = [1 -1 1 -1 1 -1]
% n=1:2:n_inf_bits;
% signal_inf_bits(n) = 1;
% 
% signal_inf_bits = ones(n_inf_bits,1);       %information signal = [-1 -1 -1 1 1 1]
% signal_inf_bits(1:fix(n_inf_bits/2)) = -1;


% adding sync marks (start)
%signal  = InsertSyncB1(signal_inf_bits,SignBarker, delay);
signal = InsertSyncB2(signal_inf_bits, sign_barker_b2, period);
signal = InsertSyncB1(signal, sign_barker_b1, delay);
% adding sync marks (stop)

samples = kt * Fs / F;       %!!!! number of samples per one symbol
if abs(samples - fix(samples)) > 0                  % check Freq assignment error
    disp(['Error. Freq assignment error. samples = kt*Fs/F = ', num2str(samples)]);
    return;
end

%x = linspace(0,kt*n_total_bits*2*pi-(F*Td),n_total_bits*samples);
x = 0:F * Td:(kt * n_total_bits * 2 * pi) - (F * Td);

signal_long = (Short2Long(signal, samples) + 1) / 2;     % OOK
%signal_long = (Short2Long(signal, samples));            % BPSK
signal_long(1:delay * samples) = 0;
signal_long = SignalLongFilter(signal_long, samples, Fs);     % filtering
u = signal_long.*sin(x)';

%signal(1:5)
% modulation(stop)

%air channel modeling (start)
u = u / std(u);
% SNR = 1000;
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

%return

%**************************************************
%***************Receiver***************************
%**************************************************

%info = audiodevinfo
%support = audiodevinfo(1,1,44100,16,1)
%Fs = 44100;% 44100; %96000 
% F = 5*Fs/100;  %frequency of signal, 200<F<Fs/2, [Hz]. even(1*Fs/100, 2*Fs/100, 4*Fs/100). F = 2(and 4)*Fs/100 - optimum
% kF = 4;         %f1 = kF*f0, kF=2(and 8) - optimum
% n_inf_bits = 1*8*1024;   %number of information bits
% SignBarker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1             1 1 1 1 1 -1 -1 1 1 -1 1 -1 1              1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
% %SignBarker = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
%tt = 1+kt*(2*length(SignBarker) + n_inf_bits)/F;   %common transmit time
n = fix(n_inf_bits / period);
tt = 1 + kt * (2 * length(sign_barker_b1) + n * length(sign_barker_b2) + n_inf_bits)/F;   %common transmit time

nBits = 24;
samples = kt * Fs / F;       %!!!! number of samples per one symbol
if abs(samples - fix(samples)) > 0                  %check Freq assignment error
    disp(['Error. Freq assignment error. samples = kt*Fs/F = ', num2str(samples)]);
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

sign_barker_b1_long = Short2Long(sign_barker_b1, samples);
sign_barker_b2_long = Short2Long(sign_barker_b2, samples);

plot_time(z, Fs, 'sec', 'recorded signal z')
plot_psd(z, Fs, 'Hz', 'PSD of received signal z');

figure, spectrogram(z, 400, 100, [], Fs); % Compute the short-time Fourier transform. Divide the waveform into 400-sample segments with 100-sample overlap
title('Received signal spectrogram');

[est_signal_b, signal_constel, index_a, index_b] = calc_ook_receiver_new(z, samples, F, Fs, sign_barker_b1_long, sign_barker_b2_long, n_inf_bits, period, signal_inf_bits);

% equalization start()
sign_x = SignalLongFilter(sign_barker_b1_long, samples, Fs);     %filtering
%sign_x = sign_barker_b1_long;
x = 0:F * Td:(kt * n_total_bits * 2 * pi) - (F * Td);
sign_x = sign_x.*sin(x(1:length(sign_x)))';
z_new = equalizer_first(sign_x, z, 3 * n_sign_barker_b1, index_a);
plot_psd(z_new, Fs, 'Hz', 'PSD of equalized received z');
% equalization stop()

[est_signal_b, signal_constel, index_a, index_b] = calc_ook_receiver_new(z_new, samples, F, Fs, sign_barker_b1_long, sign_barker_b2_long, n_inf_bits, period, signal_inf_bits);
calc_snr(z, length(sign_barker_b1_long), index_a, index_b);

%write file start
%[errmsg] = signal2file('output\output.txt', est_signal_b);
[errmsg] = signal2file(strcat('output\', file_name), est_signal_b);
%[errmsg] = signal2file(strcat('output\',char(datetime),'.txt'), est_signal_b);
if length(errmsg) ~= 0
    disp('signal2file error');
    disp(errmsg);
    return;
end
%write file stop
toc