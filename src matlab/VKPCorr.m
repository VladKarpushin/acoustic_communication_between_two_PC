%������� ������� ������� ����������� ���������� ��������
%2010.01.09 ���������� ������ � �������������� �������, ������ ����������
%������ � ����� ����� 1
function [corr] = VKPCorr(SignA,SignB)
SignA = SignA - mean(SignA);
SignB = SignB - mean(SignB);
varA = var(SignA);
varB = var(SignB);
corr = SignA.*SignB;
corr = sum(corr);
corr = corr/(sqrt(varA*varB)*(length(SignA)-1));