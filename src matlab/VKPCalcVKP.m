%функция проводит сегментацию изображения и вычисление сигнала рельефа
%добавлена проверка на длину сигнала 2010.01.08
function [Dest Err] = VKPCalcVKP(SignA,SignB,NkrVKP)
% input:
% 	SignA - указатель на первый сигнал
% 	SignB - указатель на второй сигнал
%   NkrVKP - края ВКФ
% output:
% 	Dest - указатель на сигнал, который будет хранить ВКП
%   Err - информация об ошибках
% 	Сдвигается сигнал меньшей размерности.

    Dest = 0;
    Err = false;


    %проверка на длину сигнала
    if (length(SignA) < NkrVKP) | (length(SignB) < NkrVKP) 
        '(length(SignA) < NkrVKP) | (length(SignB) < NkrVKP)'
        Err = true;
        [length(SignA) length(SignB)  NkrVKP]
    end

% определяем сигнал меньшей размерности (н)	
	if length(SignA) >= length(SignB)
		Sign1 = SignA;
		Sign2 = SignB;	%сигнал наименьшей размерности
	else
		Sign1 = SignB;
		Sign2 = SignA;
    end
% определяем сигнал меньшей размерности (к)

	BufSize1 = length(Sign1);
	BufSize2 = length(Sign2);
	BufSizeVKP = BufSize1 + BufSize2 - 1;
    Dest = zeros(1, BufSizeVKP);

% 	1-й этап
% 	до:
% 	Sign1		**************
% 	Sign2	 ****
% 	после:
% 	Sign1		**************
% 	Sign2		****
	for i = NkrVKP:BufSize2
        sign1_t = Sign1(1:i);
        sign2_t = Sign2(BufSize2-i+1:BufSize2);
		Dest(1,i) = VKPCorr(sign1_t,sign2_t);
    end

% 	2-й этап
% 	до:
% 	Sign1		**************
% 	Sign2		 ****
% 	после:
% 	Sign1		**************
% 	Sign2				  ****
    for i = BufSize2+1:BufSize1
        sign1_t = Sign1(i-BufSize2+1:i);
        sign2_t = Sign2(1:BufSize2);
		Dest(1,i) = VKPCorr(sign1_t,sign2_t);
    end

    
% 	3-й этап
% 	до:
% 	Sign1		**************
% 	Sign2				   ****
% 	после:
% 	Sign1		**************
% 	Sign2				      ****
	for i = BufSize1+1:BufSizeVKP-NkrVKP
        sign1_t = Sign1(i-BufSize2+1:BufSize1);
        sign2_t = Sign2(1:BufSize1+BufSize2-i);
		Dest(1,i) = VKPCorr(sign1_t,sign2_t);
    end
%     length(sign1_t)