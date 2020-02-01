% Function calculates cross correlation function (CCF)

function [ccf err] = calc_ccf_fft(sign_a, sign_b, n_border)
% input:
% 	sign_a - ��������� �� ������ ������
% 	sign_b - ��������� �� ������ ������ (��������, �������)
% output:
% 	ccf - ��������� �� ������, ������� ����� ������� ���
%   err - ���������� �� �������
% 	���������� ������ ������� �����������.

ccf = 0;
err = false;
%�������� �� ����� �������
if (length(sign_a) < n_border) | (length(sign_b) < n_border) 
    disp('(length(sign_a) < n_border) | (length(sign_b) < n_border)');
    err = true;
    [length(sign_a) length(sign_b)  n_border]
end

% ���������� ������ ������� ����������� (�)	
% sign1 = sign_b;
% sign2 = sign_a; % ������ ���������� �����������
% if length(sign_a) >= length(sign_b)
%     sign1 = sign_a;
%     sign2 = sign_b;	% ������ ���������� �����������
% end

if length(sign_b) > length(sign_a)
    tmp = sign_a;
    sign_a = sign_b;    % long signal
    sign_b = tmp;       % short signal
end
% ���������� ������ ������� ����������� (�)

% ������ ������� ������ ����������� (�)
% tmp = sign2;
% sign2 = zeros(length(sign1), 1);
% sign2(1:length(tmp)) = tmp;
tmp = sign_b;
%sign_b = zeros(length(sign_a), 1);
sign_b = zeros(size(sign_a));
sign_b(1:length(tmp)) = tmp;

% ������ ������� ������ ����������� (�)

sign_a_fft = fft(sign_a);
sign_b_fft = fft(sign_b);
ccf = real(ifft(sign_a_fft .* conj(sign_b_fft)));
den = max(sign_a) * max(sign_b) * length(tmp); % denominator
ccf = ccf / den;
end