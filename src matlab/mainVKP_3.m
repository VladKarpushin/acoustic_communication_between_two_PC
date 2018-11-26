%2016-02-15
%Tool for synchronization sequence choosing
close all,clc,clear all;

SignA = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1         -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1       1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]';    %Barker code
%SignA = [1 -1 -1 -1 -1 1 -1 -1 -1 1 -1 -1 -1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 1]';  % Gold
SignB = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code
T=100;

%SignB = SignA; 
[Dest Err] = VKPCalcVKP(SignA,SignB,2);
figure,plot(Dest);
title('cross-correlation function (SignA,SignB)');

SignALong = ones(length(SignA)*T,1);
SignBLong = ones(length(SignB)*T,1);
for n = 1:length(SignA)
    iA = 1 + (n-1)*T;
    iB = iA+T;
    if SignA(n) == -1
        SignALong(iA:iB) = -1;
    end
end
for n = 1:length(SignB)
    iA = 1 + (n-1)*T;
    iB = iA+T;
    if SignB(n) == -1
        SignBLong(iA:iB) = -1;
    end
end

[DestLong Err] = VKPCalcVKP(SignALong,SignBLong,2);

[Max,Imax] = max(Dest)              %largest element index
[MaxLong,ImaxLong] = max(DestLong)  %largest element index

figure,plot(DestLong);
title('cross-correlation function (SignALong,SignBLong)');

figure,stem(SignA);
title('SignA');

figure,stem(SignALong);
title('SignALong');

figure,stem(SignB);
title('SignB');

figure,stem(SignBLong);
title('SignBLong');