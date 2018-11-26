close all,clc,clear all;

%SignA = [1 -1 1 -1 1 -1 1 -1 1 -1 ];
SignA = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1         -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1       1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];    %Barker code
%SignA = [1    -1    -1     1    -1     1     1    -1     1    -1    -1 1 -1     1     1    -1]; %hadamard
%SignA = [1 0 0 0 0 1 0 0 0 1 0 0 0 1 0 1 0 0 0 1 1 0 0 0 1 1 0 1 0 1 1];  % Gold
SignB = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]; %Barker code
%SignB = SignA; 
[Dest Err] = VKPCalcVKP(SignA,SignB,2);
% x = 1:floor(length(Dest));
% x=x/Fs;
% figure,plot(x,Dest);
% xlabel('sec');
figure,plot(Dest);
title('cross-correlation function');

T=10;
SignALong = zeros(length(SignA)*T,1);
for n = 1:length(SignA)
    SignALong()
    iA = 1 + (n-1)*Samples;
    iB = iA+Samples;
    SignALong(iA:iB) = 1;
end
