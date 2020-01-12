% 2019-12-08
% signal-to-noise (snr) calculation
% 2020-01-02
% took into account pll block

function [snr_estimated] = calc_snr(z, ind_a, ind_b, freq_burst_size)

index_noise = ind_a - 7 - freq_burst_size;
if (index_noise > 1) && (ind_b > 1) && (ind_a < length(z)) && (ind_b <= length(z))  %snr estimation 
    s = std(z(ind_a:ind_b));
    n = std(z(1:index_noise));
    snr_estimated = (s / n) - 1; % maybe "snr_estimated = (s / n);" is correct
    disp(['snr estimated = ', num2str(round(snr_estimated)), ' [times]']);
    disp(['snr estimated = ', num2str(round(10 * log10(snr_estimated))), ' [dB]']);
else
    disp('snr estimation error');
end