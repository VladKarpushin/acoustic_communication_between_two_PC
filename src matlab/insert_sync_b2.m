% function [Dest] = InsertSyncB2(Signlnf,SignSync, period)
% this function inserts SignSyncB2 into Signlnf with period

function [sign_out] = insert_sync_b2(sign_inf,sign_sync, period)
% input:
% sign_inf   - information signal
% sign_sync  - sync signal B2
% period    - period of B2 insertation
% output:
% sign_out      - information signal with B2 sync signal every period bits

n = fix(length(sign_inf) / period);
sign_out = zeros(length(sign_inf) + length(sign_sync) * n, 1) ;
for i = 1:n
    iA = 1 + (i - 1) * period;
    iB = iA + period - 1;
    iAA = 1 +	(i - 1) * (period + length(sign_sync));
    iBB = iAA + period - 1;
    sign_out(iAA:iBB) = sign_inf(iA:iB);
    sign_out(iBB+1:iBB + length(sign_sync)) = sign_sync;
end
if abs((length(sign_inf) / period) - n) > 0
    iA = 1 + n * period;
    iB = length(sign_inf) ;
    iAA = 1 + n * (period + length(sign_sync));
    iBB = length(sign_out);
    sign_out(iAA:iBB) = sign_inf(iA:iB);
end
end