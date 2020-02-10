%function [est_signal_b Err] = CalcSignalEstimation(CorrIntegrall,threshold)
%This function estimates information bits (information signal)
%without plots
%for BPSK only
%2016-11-05 added code for signal constellation
%2016-12-20 added iAB1 and iBB1 export for SNR calculation
function [est_signal_b, max_sync_b1, max_sync_b2, err, signal_constel, ind_a, ind_b] = calc_signal_estimation_bpsk(threshold, sign_barker_b1_long, sign_barker_b2_long, samples, signal_complex)
% input:
% 	corr_integral    - cross correlation function (CCF) received signal SignR and sin wave. Another name is correlation integral
% 	threshold       - resolver threshold. Should be zero for BPSK
%   sign_barker_long  -    sync signal(long)
%   samples         - quantity of samples per one symbol
%   signal_complex   - complex signal
% output:
% 	est_signal_b     - estimated information bits
%   err             - error information
%   max_sync_b1     - max_sync_b1 is max(CCF)
%   max_sync_b2     - max_sync_b2 is min(CCF)
%   delta           - delta is difference between index(max_sync_b2) and index(max_sync_b1)
%   StdSignSync     - StdSignSync is std(CCF) between two mainlobes
%   signal_constel   - signal constellation
%   ind_a            - index of first symbol. For SNR estimation
%   ind_b            - index of last symbol. For SNR estimation

err = 0;
max_sync_b1 = 0;
max_sync_b2 = 0;
est_signal_b = 0;
signal_constel = 0;
ind_a = 0;
ind_b = 0;

%corr_integral = real(signal_complex);
est_signal_long = zeros(length(signal_complex), 1);
est_signal_long = (2 * (real(signal_complex) > threshold)) - 1;    %resolver
est_signal_long_q = (2 * (imag(signal_complex) > threshold)) - 1;    % resolver for q part of qpsk

%****syncronization start*******
sync_b1 = calc_ccf_fft(est_signal_long, sign_barker_b1_long, 0);
sync_b2 = calc_ccf_fft(est_signal_long, sign_barker_b2_long, 0);

[max_sync_b1, ind_max_sync_b1] = max(abs(sync_b1));
[max_sync_b2, ind_max_sync_b2] = max(abs(sync_b2));

%if abs(min(sync_b1)) > abs(max(sync_b1))
if sign(sync_b1(ind_max_sync_b1)) < 0
    est_signal_long = -est_signal_long;
end

ind_a = ind_max_sync_b1 + length(sign_barker_b1_long);
ind_b = ind_max_sync_b2 - 1;

est_signal_b = long_to_short(est_signal_long(ind_a:ind_b), samples);
est_signal_b_q = long_to_short(est_signal_long_q(ind_a:ind_b), samples);
signal_constel = long_to_short(signal_complex(ind_a:ind_b), samples);
%****syncronization stop*******


if abs(length(est_signal_b) / 8 - fix(length(est_signal_b) / 8)) > 0                  %checking if length of est_signal_b is multiple with 8
    disp(['Error. abs(length(est_signal_b)/8 - fix(length(est_signal_b)/8)) > 0. length(est_signal_b) = ', num2str(length(est_signal_b))]);
    err = 1;
    return;
end
end
