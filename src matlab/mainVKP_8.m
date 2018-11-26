%2016-02-15
%Tool for synchronization sequence choosing
close all,clc,clear all;

%SignA = [ 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1         -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1       1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]';    %Barker code
%SignA = [1 -1 -1 -1 -1 1 -1 -1 -1 1 -1 -1 -1 1 -1 1 -1 -1 -1 1 1 -1 -1 -1 1 1 -1 1 -1 1 1]';  % Gold

%Barker code N=13
%SignA = [ 1 -1 -1 1         1 1 1 1 1 -1 -1 1 1 -1 1 -1 1         -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1       1 1 1 1 1 -1 -1 1 1 -1 1 -1 1         1 1 1 0 0 1 0 1]';    %Barker code N=13
%SignB = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
% SignA = [ 1 1 1 1 1 1       1 1 1 1 1 -1 -1 1 1 -1 1 -1 1         1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1     1 1 1 1 1 -1 -1 1 1 -1 1 -1 1    1 1 1 -1 -1 1 -1 1]';    %Barker code N=13
% SignB = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1    1 1 1 1 1 -1 -1 1 1 -1 1 -1 1          1 1 1 1 1 -1 -1 1 1 -1 1 -1 1        1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
SignA = [ 1 1 1 1 1 1       1 1 1 1 1 -1 -1 1 1 -1 1 -1 1   1 1 1 -1 -1 1 -1 1      -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1]';    %Barker code N=13
SignB = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13

%Barker code N=13 (4 times)
% SignA = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1            1 1 1 1 1 -1 -1 1 1 -1 1 -1 1        1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13
% SignB = SignA;

%Barker code N=11
% SignA = [1 -1 -1 1 1 -1      1 1 1 -1 -1 -1 1 -1 -1 1 -1       -1 -1 -1 1 1 1 -1 1 1 -1 1      1 1 1 -1 -1 -1 1 -1 -1 1 -1          -1 -1 -1 1 1 1 -1 1 1 -1 1          -1 -1 -1 1 1 1 -1 1 1 -1 1                        1 1 -1 1 1 -1 -1]';    %Barker code N=11
% SignB = [1 1 1 -1 -1 -1 1 -1 -1 1 -1]'; %Barker code N=11

T=100;

%SignB = SignA; 
[Dest Err] = VKPCalcVKP(SignA,SignB,length(SignB));
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

tic
[DestLong Err] = VKPCalcVKP(SignALong,SignBLong,length(SignBLong));
toc

[DestFFT Err] = VKPCalcVKP_FFT(SignA,SignB,0);

tic
[DestLongFFT Err] = VKPCalcVKP_FFT(SignALong,SignBLong,0);
toc

[Max,Imax] = max(DestFFT)              %largest element index
[MaxLong,ImaxLong] = max(DestLongFFT)  %largest element index

[Min,Imin] = min(DestFFT)  %smalles element index
[MinLong,IminLong] = min(DestLongFFT)  %smalles element index

iA = Imax + length(SignB);
iB = Imin-1;
SignA(iA:iB)

iALong = ImaxLong + length(SignBLong);
iBLong = IminLong-1;
tmp = SignALong(iALong:iBLong);


figure,plot(DestFFT);
title('cross-correlation function DestFFT');

figure,plot(DestLongFFT);
title('cross-correlation function DestLongFFT');

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