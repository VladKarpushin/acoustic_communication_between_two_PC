% function converts signal from short form to long. For example, short
% singnal is [-1, 1], samples = 2, then long signal is [-1 -1, 1  1]
% this is interpolator

function [sign_long] = short_to_long(sign_short, Samples)
% input:
% sign_short - short signal
% Samples
% output:
% sign_long - long signal
sign_long = zeros(length(sign_short) * Samples, 1);
for n = 1:length(sign_short)
    iA = 1 + (n - 1) * Samples;
    iB = iA + Samples - 1;
    if sign_short(n) == -1
        sign_long(iA:iB) = -1;
    elseif sign_short(n) == 1
        sign_long(iA:iB) = 1;        
    end
end