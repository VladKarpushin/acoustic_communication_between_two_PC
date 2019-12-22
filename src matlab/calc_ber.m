% 2019-12-22
% bit error rate calculation (BER)

function BER = calc_ber(sing_orig, sing_est, n_inf_bits)

BER = -2;
if length(sing_est) == n_inf_bits
    BER = mean(abs(sing_orig - sing_est) / 2);   %The bit error rate (BER) calculation
else
    disp(['Error. Can not calculate BER, because of different array size. length(est_signal_b) = ', num2str(length(sing_est)), ',  length(signal_inf_bits) = ', num2str(length(sing_orig))]);
end