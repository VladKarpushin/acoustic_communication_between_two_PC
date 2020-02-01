% function converts signal from short form to long. For example, short
% singnal is [-1, 1], samples = 2, then long signal is [-1 -1, 1  1]
% this is interpolator function

function [sign_long] = short_to_long(sign_short, samples)
% input:
% sign_short - short signal
% samples
% output:
% sign_long - long signal
sign_long = zeros(length(sign_short) * samples, 1);
for n = 1:length(sign_short)
    ind_a = 1 + (n - 1) * samples;
    ind_b = ind_a + samples - 1;
    if sign_short(n) == -1
        sign_long(ind_a:ind_b) = -1;
    elseif sign_short(n) == 1
        sign_long(ind_a:ind_b) = 1;        
    end
end