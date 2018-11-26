%this function inserts SignSyncB1 into beginning and end of Signlnf
function [Dest] = InsertSyncB1(Signlnf,SignSync, delay)
% input:
% Signlnf   - information signal
% SignSync  - sync signal B1
% delay     - time delay in a beginning of transmission (unit is bit)
% output:
% Dest      - information signal with B1 sync signal in beginning and end
% of signal


nTotalBits = delay + 2*length(SignSync) + length(Signlnf);
Dest = zeros(nTotalBits, 1) ;
Dest(1+delay:length(SignSync)+delay) = SignSync;
Dest(nTotalBits - length(SignSync)+1:nTotalBits) = -SignSync;
Dest(length(SignSync) + 1 + delay:length(SignSync) + delay + length(Signlnf)) = Signlnf;
% 
% signal(1+delay:length(SignBarker)+delay) = SignBarker;
% signal(nTotalBits - length(SignBarker)+1:length(signal)) = -SignBarker;
% signal(length(SignBarker)+1+delay:length(SignBarker)+delay + nInfBits) = signalInf_b;

end