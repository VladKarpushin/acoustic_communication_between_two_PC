%function converts signal from long form to short. For example, long signal is [-1 -1, 1  1], samples = 2, then short singnal is [-1, 1] 
% This is sampler function

function [sign_short] = long_to_short(sign_long, samples)
% input:
% sign_long - long signal
% samples
% output:
% SignShort - short signal

% sign_long = ones(length(SignShort)*samples,1);
% for n = 1:length(SignShort)
%     iA = 1 + (n-1)*samples;
%     iB = iA+samples-1;
%     if SignShort(n) == -1
%         sign_long(iA:iB) = -1;
%     end
% end


ind_a = fix(samples / 2);
ind_b = length(sign_long);
short_size = round(ind_b / samples);
sign_short =  zeros(short_size, 1);
i = 1;
for n = ind_a:samples:ind_b
    sign_short(i) = sign_long(n);
    i = i + 1;
end