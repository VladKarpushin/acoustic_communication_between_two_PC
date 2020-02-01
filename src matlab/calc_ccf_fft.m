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
sign1 = sign_b;
sign2 = sign_a; % ������ ���������� �����������
if length(sign_a) >= length(sign_b)
    sign1 = sign_a;
    sign2 = sign_b;	% ������ ���������� �����������
end
% ���������� ������ ������� ����������� (�)

% ������ ������� ������ ����������� (�)
sign_tmp = sign2;
sign2 = zeros(length(sign1), 1);
sign2(1:length(sign_tmp)) = sign_tmp;
% ������ ������� ������ ����������� (�)

sign1_fft = fft(sign1);
sign2_fft = fft(sign2);
ccf = real(ifft(sign1_fft .* conj(sign2_fft)));
den = max(sign_a) * max(sign_b) * length(sign_b);
ccf = ccf / den;
end