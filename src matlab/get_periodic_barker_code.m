% Function returns periodic Barker code

function [out] = get_periodic_barker_code(barker_code, n_period)
% input:
% 	barker_code - single Barker code
% 	n_period - quantity of periods
% output:
% 	out - periodic Barker code

out = zeros(length(barker_code) * n_period, 1);
for n = 1:n_period
    ind_a = 1 + (n - 1) * length(barker_code);
    ind_b = ind_a + length(barker_code) - 1;
    out(ind_a:ind_b) = barker_code;
end
end
