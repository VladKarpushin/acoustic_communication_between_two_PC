%2016-07-27
%This script plots level of CCF sidelobe  from length of sync sequence
% x - quantity of Barker codes. x*13 - length of sync sequence
% y - level of CCF sidelobe. It was found by std function
close all,clc,clear all;
y = [0.2330    0.1650    0.1350    0.1180    0.0990    0.0930    0.0860    0.0820    0.0730    0.0710    0.0730    0.0750    0.0630    0.0570    0.0500    0.0540    0.0510    0.0580    0.0510 0.0500    0.0400    0.0310    0.0380    0.0380    0.0340    0.0240    0.0290    0.0200    0.0200    0.0270    0.0160];
x = [1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18    19    20    25    30    35    40    45    50    55    60    65    70    75];

figure, plot(x,y,'o',x,1./(x*13),'r',x,1./sqrt(x*13),'b');
%figure, plot(x,y,'o',x,1./sqrt(x*13));
title('StdSignSync (level of CCF sidelobe)');
xlabel('nSignBarker (quantity of Barker codes N=13)');
ylim([0 0.25]);