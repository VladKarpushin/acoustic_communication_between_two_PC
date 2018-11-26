%function converts signal from long form to short. For example, long signal is [-1 -1, 1  1], samples = 2, then short singnal is [-1, 1] 

function [SignShort] = Long2Short(SignLong, Samples)
% input:
% SignLong - long signal
% Samples
% output:
% SignShort - short signal

% SignLong = ones(length(SignShort)*Samples,1);
% for n = 1:length(SignShort)
%     iA = 1 + (n-1)*Samples;
%     iB = iA+Samples-1;
%     if SignShort(n) == -1
%         SignLong(iA:iB) = -1;
%     end
% end


iA = fix(Samples/2);
iB = length(SignLong);
nMax = round(iB/Samples);
SignShort =  zeros(nMax,1);
i = 1;
for n = iA:Samples:iB
    SignShort(i) = SignLong(n);
    i = i + 1;
end
a=1;