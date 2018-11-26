%2016-08-06
%sync test with Barker codes of different length

close all,clc,clear all;


%*****Barker codes set generation (start)*****
nSignBarker = 75;   %quantity of Barker codes in a set.
%SignBarkerOneB1 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems. Barker codes have length at most 13 and have low correlation sidelobes
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
title('ccfB1B1');

ccfB2B2 = CalcCCF_FFT(SignBarkerB2,SignBarkerB2, 0);
figure, plot(ccfB2B2);
title('ccfB2B2');

ccfB1B2 = CalcCCF_FFT(SignBarkerB1,SignBarkerB2, 0);
figure, plot(ccfB1B2);
title('ccfB1B2');