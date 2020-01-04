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
%2016-08-21 added sign_barker_long = Short2Long(sign_barker, samples);
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
[signal_inf_bits, errmsg] = file2signal(strcat('..\input\', filename));

if ~isempty(errmsg)
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
F = Fs / 30;  %frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave
%F = Fs/5;  %frequency of signal, 200<F<Fs/2, [Hz]. F = Fs/14 - max, F = Fs/30 - max for Fs = 96000; For example, F = Fs/30, 30 - number of samples per one wave
kt = 2;     %coefficient of duration of one symbol, kt/F = duration of one symbol
n_inf_bits = 1024 * 4 * 1;      % number of information bits
%n_inf_bits = length(signal_inf_bits);
Td = 2 * pi / Fs;   % sampling interval
delay = 1000;       % time delay in a beginning of transmission (unit is bit)
pll_block_size = 1000;

%*****Barker codes set generation (start)*****
n_sign_barker = 75;   %quantity of Barker codes in a set.
sign_barker_one = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems. Barker codes have length at most 13 and have low correlation sidelobes
sign_barker = GetPeriodicBarkerCode(sign_barker_one, n_sign_barker);
%*****Barker codes set generation (stop)*****

n_total_bits = delay + pll_block_size + 2 * length(sign_barker) + n_inf_bits;
show_sign_para(kt, F, Fs, n_total_bits, n_inf_bits, delay);
%return

%modulation(start)
signal_inf_bits = 2 * randi([0, 1], n_inf_bits, 1) - 1; % model of information signal is noise

% signal_inf_bits = -1*ones(n_inf_bits,1);       %information signal = [1 -1 1 -1 1 -1]
% n=1:2:n_inf_bits;
% signal_inf_bits(n) = 1;
% 
% signal_inf_bits = ones(n_inf_bits,1);       %information signal = [-1 -1 -1 1 1 1]
% signal_inf_bits(1:fix(n_inf_bits/2)) = -1;

%adding sync marks (start)
%signal  = InsertSyncB1(signal_inf_bits, sign_barker, delay);
signal  = construct_signal_bpsk(signal_inf_bits, sign_barker, delay, pll_block_size);
%adding sync marks (stop)

samples = kt * Fs / F;       %!!!! number of samples per one symbol
if abs(samples - fix(samples)) > 0                  % check Freq assignment error
    disp(['Error. Freq assignment error. samples = kt*Fs/F = ', num2str(samples)]);
    return;
end

x = 0:F * Td:(kt * n_total_bits * 2 * pi) - (F * Td);
%x = linspace(0,kt*n_total_bits*2*pi-(F*Td),n_total_bits*samples);
%signal_long = (Short2Long(signal, samples)+1)/2;     %OOK
signal_long = Short2Long(signal, samples);        %BPSK
%signal_long(1:delay * samples) = 0;
signal_long = SignalLongFilter(signal_long, samples, Fs);     %filtering

% qpsk_part = 2 * randi([0, 1], length(signal), 1) - 1; % model of information signal is noise
% qpsk_part = (Short2Long(qpsk_part, samples));        %BPSK
% qpsk_part = SignalLongFilter(qpsk_part, samples, Fs);     %filtering
%u = signal_long.*sin(x)' + qpsk_part.*cos(x)';

u = signal_long .* sin(x)';
%modulation(stop)

%air channel modeling (start)
u = u / std(u);
%SNR = 10;        %signal to noise ratio
%signal_noise = randn(length(u), 1) / SNR;
%u = u + signal_noise;
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

tt = 1 + kt * (pll_block_size + 2 * length(sign_barker) + n_inf_bits) / F;   %common transmit time
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
recordblocking(recObj, tt + 2);
disp('End of Recording.');

% Store data in double-precision array.
z = getaudiodata(recObj)';      %received signal
%z = u';

sign_barker_long = Short2Long(sign_barker, samples);

plot_time(z, Fs, 'sec', 'recorded signal z')
plot_psd(z, Fs, 'Hz', 'PSD of received signal z');

figure, spectrogram(z, 400, 100, [], Fs); % Compute the short-time Fourier transform. Divide the waveform into 400-sample segments with 100-sample overlap
title('Received signal spectrogram');

%Fs = Fs + Fs * 7 * 10^-6; % actual drift is 0.1544 Hz
[~, index_a, ~] = calc_bpsk_receiver(z, samples, F, Fs, sign_barker_long, n_inf_bits, signal_inf_bits);

% equalizer start()
sign_x = SignalLongFilter(sign_barker_long, samples, Fs);     %filtering
%sign_x = sign_barker_long;
x = 0:F * Td:(kt * n_total_bits * 2 * pi) - (F * Td);
sign_x = sign_x .* sin(x(1:length(sign_x)))';
z_new = equalizer_first(sign_x, z, 3 * n_sign_barker, index_a);
plot_psd(z_new, Fs, 'Hz', 'PSD of equalized received z');
% equalizer stop()

[est_signal_b, index_a, index_b] = calc_bpsk_receiver(z_new, samples, F, Fs, sign_barker_long, n_inf_bits, signal_inf_bits);
calc_snr(z, length(sign_barker_long), index_a, index_b, pll_block_size * samples);
est_F = calc_freq_offset(z, length(sign_barker_long), index_a, pll_block_size * samples, Fs);
disp(['est_F - F = ', num2str(est_F - F)]);

%(F - est_F) * 7 - Fs freq offset

%write file (start)
[errmsg] = signal2file('output\output.txt', est_signal_b);
%[errmsg] = signal2file(strcat('output\',char(datetime),'.txt'), est_signal_b);
if ~isempty(errmsg)
    disp('signal2file error');
    disp(errmsg);
    return;
end
%write file (stop)