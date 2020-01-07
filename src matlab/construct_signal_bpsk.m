% this function constructs signal for bps transmission

function [sign_out] = construct_signal_bpsk(sign_inf, sign_sync, delay, freq_burst_size)
% input:
% sign_inf   - information signal
% sign_sync  - sync signal B1
% delay     - time delay in a beginning of transmission (unit is bit)
% output:
% sign_out      - information signal with B1 sync signal in beginning and end
% of signal


n_total_bits = delay + freq_burst_size + 2 * length(sign_sync) + length(sign_inf);
sign_out = zeros(n_total_bits, 1);
sign_out(1 + delay + freq_burst_size:delay + freq_burst_size + length(sign_sync)) = sign_sync;
sign_out(1 + delay:delay + freq_burst_size) = 1;
sign_out(n_total_bits - length(sign_sync) + 1:n_total_bits) = - sign_sync;
sign_out(1 + delay + freq_burst_size + length(sign_sync):delay + freq_burst_size + length(sign_sync) + length(sign_inf)) = sign_inf;
end