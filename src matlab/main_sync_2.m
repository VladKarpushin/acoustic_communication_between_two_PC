%2016-08-06
%sync test with Barker codes of different length

close all,clc,clear all;


nSignBarker = 75;   %quantity of Barker codes in a set.
SignBarkerOneB1 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems. Barker codes have length at most 13 and have low correlation sidelobes
SignBarkerOneB2 = [1 1 1 -1 -1 -1 1 -1 -1 1 -1]'; %Barker code N=11. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems.

SignBarkerB1 = GetPeriodicBarkerCode(SignBarkerOneB1, nSignBarker);
SignBarkerB2 = GetPeriodicBarkerCode(SignBarkerOneB2, nSignBarker);

ccfB1B1 = CalcCCF_FFT(SignBarkerOneB1,SignBarkerOneB1, 0);
figure, plot(ccfB1B1);
title('ccfB1B1 one');

ccfB2B2 = CalcCCF_FFT(SignBarkerOneB2,SignBarkerOneB2, 0);
figure, plot(ccfB2B2);
title('ccfB2B2 one');

ccfB1B2 = CalcCCF_FFT(SignBarkerOneB1,SignBarkerOneB2, 0);
figure, plot(ccfB1B2);
title('ccfB1B2 one');

ccfB1B1 = CalcCCF_FFT(SignBarkerB1,SignBarkerB1, 0);
figure, plot(ccfB1B1);
title('ccfB1B1 periodic');

ccfB2B2 = CalcCCF_FFT(SignBarkerB2,SignBarkerB2, 0);
figure, plot(ccfB2B2);
title('ccfB2B2 periodic');

ccfB1B2 = CalcCCF_FFT(SignBarkerB1,SignBarkerB2, 0);
figure, plot(ccfB1B2);
title('ccfB1B2 periodic');


period = 2048;
nInfBits = 2048*4;
signalInf_b = 2*randi([0,1],nInfBits,1)-1; %information signal = noise
% signalInf_b = 1:14;
% SignSyncB2      = [-2 -2];
% SignSyncB1      = [-1 -1 -1];
signalInf_B2    = InsertSyncB2(signalInf_b,3*SignBarkerB2, period);
signalInf_B2B1  = InsertSyncB1(signalInf_B2,2*SignBarkerB1);

figure, plot(signalInf_B2B1);
title('signalInf_B2B1');

ccfInfB1 = CalcCCF_FFT(signalInf_B2B1,SignBarkerB1, 0);
figure, plot(ccfInfB1);
title('ccfInfB1');

ccfInfB2 = CalcCCF_FFT(signalInf_B2B1,SignBarkerB2, 0);
figure, plot(ccfInfB2);
title('ccfInfB2');
