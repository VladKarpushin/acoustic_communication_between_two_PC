% 2019-12-08
% signal-to-noise (snr) calculation

function [snr_estimated] = calc_snr(z, len_of_barker, index_a, index_b)

index_a = index_a - len_of_barker;
index_b = index_b + len_of_barker;
if (index_a-4 > 1) && (index_b > 1) && (index_a < length(z)) && (index_b < length(z))  %snr estimation 
    s = std(z(index_a:index_b));
    n = std(z(1:index_a - 4));
    snr_estimated = (s / n) - 1;
    disp(['snr estimated = ', num2str(round(snr_estimated)), ' [times]']);
    disp(['snr estimated = ', num2str(round(10 * log10(snr_estimated))), ' [dB]']);
else
    disp('snr estimation error');
end