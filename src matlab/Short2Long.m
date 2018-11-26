%function converts signal from short form to long. For example, short
%singnal is [-1, 1], samples = 2, then long signal is [-1 -1, 1  1]

function [SignLong] = Short2Long(SignShort, Samples)
% input:
% SignShort - short signal
% Samples
% output:
% SignLong - long signal
SignLong = ones(length(SignShort)*Samples,1);
for n = 1:length(SignShort)
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples-1;
    if SignShort(n) == -1
        SignLong(iA:iB) = -1;
    end
end