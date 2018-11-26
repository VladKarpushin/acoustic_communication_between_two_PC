%функция функция находит коэффициент корреляции векторов
%2010.01.09 исправлена ошибка в корреляционной функции, теперь корреляция
%самого с собой равна 1
function [corr] = VKPCorr(SignA,SignB)
SignA = SignA - mean(SignA);
SignB = SignB - mean(SignB);
varA = var(SignA);
varB = var(SignB);
corr = SignA.*SignB;
corr = sum(corr);
corr = corr/(sqrt(varA*varB)*(length(SignA)-1));