% this function inserts sign_syncB1 into a beginning and an end of signal
function [sign_out] = InsertSyncB1(sign_inf, sign_sync, delay)
% input:
% sign_inf   - information signal
% sign_sync  - sync signal B1
% delay     - time delay in a beginning of transmission (unit is bit)
% output:
% sign_out      - information signal with B1 sync signal in beginning and end
% of signal


n_total_bits = delay + 2*length(sign_sync) + length(sign_inf);
sign_out = zeros(n_total_bits, 1);
sign_out(1 + delay:length(sign_sync) + delay) = sign_sync;
sign_out(n_total_bits - length(sign_sync) + 1:n_total_bits) = - sign_sync;
sign_out(length(sign_sync) + 1 + delay:length(sign_sync) + delay + length(sign_inf)) = sign_inf;
end