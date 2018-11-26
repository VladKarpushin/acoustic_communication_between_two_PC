%Function calculates cross correlation function (CCF)
function [Dest Err] = CalcCCF_FFT(SignA,SignB, NkrVKP)
% input:
% 	SignA - ��������� �� ������ ������
% 	SignB - ��������� �� ������ ������ (��������, �������)
% output:
% 	Dest - ��������� �� ������, ������� ����� ������� ���
%   Err - ���������� �� �������
% 	���������� ������ ������� �����������.

Dest = 0;
Err = false;
%�������� �� ����� �������
if (length(SignA) < NkrVKP) | (length(SignB) < NkrVKP) 
    '(length(SignA) < NkrVKP) | (length(SignB) < NkrVKP)'
    Err = true;
    [length(SignA) length(SignB)  NkrVKP]
end

% ���������� ������ ������� ����������� (�)	
if length(SignA) >= length(SignB)
    Sign1 = SignA;
    Sign2 = SignB;	%������ ���������� �����������
else
    Sign1 = SignB;
    Sign2 = SignA;
end
% ���������� ������ ������� ����������� (�)

Sign_tmp = Sign2;
Sign2 = zeros(length(Sign1),1);
Sign2(1:length(Sign_tmp)) = Sign_tmp;

Sign1FFT = fft(Sign1);
Sign2FFT = fft(Sign2);
Dest = real(ifft(Sign1FFT.*conj(Sign2FFT)));
Den = max(SignA)*max(SignB)*length(SignB);
Dest = Dest/Den;
end

