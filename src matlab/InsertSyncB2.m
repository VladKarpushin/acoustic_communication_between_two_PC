%function [Dest] = InsertSyncB2(Signlnf,SignSync, period)
%this function inserts SignSyncB2 into Signlnf with period
function [Dest] = InsertSyncB2(Signlnf,SignSync, period)
% input:
% Signlnf   - information signal
% SignSync  - sync signal B2
% period    - period of B2 insertation
% output:
% Dest      - information signal with B2 sync signal every period bits

n = fix(length(Signlnf)/period);
Dest = zeros(length(Signlnf) + length(SignSync)*n, 1) ;
for i = 1:n
    iA = 1 + (i-1)*period;
    iB = iA + period-1;
    iAA = 1 +	(i-1)*(period + length(SignSync)) ;
    iBB = iAA + period-1;
    Dest(iAA:iBB) = Signlnf(iA:iB);
    Dest(iBB+1:iBB+length(SignSync)) = SignSync;
end
if abs((length(Signlnf)/period) - n) > 0
    iA = 1 + n*period;
    iB = length(Signlnf) ;
    iAA = 1 + n*(period + length(SignSync));
    iBB = length(Dest);
    Dest(iAA:iBB) = Signlnf(iA:iB) ;
end

end