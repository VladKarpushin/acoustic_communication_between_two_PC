% 2019-12-22
% bit error rate calculation (BER)

function BER = calc_ber(sign_orig, sign_est, n_inf_bits, m_type)

if nargin < 4
    m_type = 'normal'; 
end

BER = -2;

if  m_type == 'debug'
    BER = mean(abs(sign_orig(1:n_inf_bits) - sign_est(1:n_inf_bits)) / 2);   %The bit error rate (BER) calculation
end

if length(sign_est) == n_inf_bits
    BER = mean(abs(sign_orig - sign_est) / 2);   %The bit error rate (BER) calculation
else
    disp(['Error. Can not calculate BER, because of different array size. length(est_signal_b) = ', num2str(length(sign_est)), ',  length(signal_inf_bits) = ', num2str(length(sign_orig))]);
end