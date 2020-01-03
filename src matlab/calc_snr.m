% 2019-12-08
% signal-to-noise (snr) calculation
% 2020-01-02
% took into account pll block

function [snr_estimated] = calc_snr(z, len_of_barker, index_a, index_b, pll_block_size)

index_a = index_a - len_of_barker;
index_b = index_b + len_of_barker;
index_noise = index_a - 7 - pll_block_size;
if (index_noise > 1) && (index_b > 1) && (index_a < length(z)) && (index_b <= length(z))  %snr estimation 
    s = std(z(index_a:index_b));
    n = std(z(1:index_noise));
    snr_estimated = (s / n) - 1; % maybe "snr_estimated = (s / n);" is correct
    disp(['snr estimated = ', num2str(round(snr_estimated)), ' [times]']);
    disp(['snr estimated = ', num2str(round(10 * log10(snr_estimated))), ' [dB]']);
else
    disp('snr estimation error');
end