%Function returns periodic Barker code
function [Dest] = GetPeriodicBarkerCode(BarkerCode,nPeriod)
% input:
% 	BarkerCode - single Barker code
% 	nPeriod - quantity of periods
% output:
% 	Dest - periodic Barker code

Dest = zeros(length(BarkerCode)*nPeriod,1);
for n = 1:nPeriod
    iA = 1 + (n-1)*length(BarkerCode);
    iB = iA + length(BarkerCode)-1;
    Dest(iA:iB) = BarkerCode;
end
end
