%2016-08-06
%sync test with Barker codes of different length
%2016-08-13. added InsertSyncB2 and InsertSyncB1
%2016-08-20 recoverred inf signal

close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

nSignBarker = 75;   %quantity of Barker codes in a set.
SignBarkerOneB1 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems. Barker codes have length at most 13 and have low correlation sidelobes
SignBarkerOneB2 = [1 1 1 -1 -1 -1 1 -1 -1 1 -1]'; %Barker code N=11. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems.

SignBarkerB1 = GetPeriodicBarkerCode(SignBarkerOneB1, nSignBarker);
SignBarkerB2 = GetPeriodicBarkerCode(SignBarkerOneB2, nSignBarker);

% ccfB1B1 = CalcCCF_FFT(SignBarkerOneB1,SignBarkerOneB1, 0);
% figure, plot(ccfB1B1);
% title('ccfB1B1 one');
% 
% ccfB2B2 = CalcCCF_FFT(SignBarkerOneB2,SignBarkerOneB2, 0);
% figure, plot(ccfB2B2);
% title('ccfB2B2 one');
% 
% ccfB1B2 = CalcCCF_FFT(SignBarkerOneB1,SignBarkerOneB2, 0);
% figure, plot(ccfB1B2);
% title('ccfB1B2 one');
% 
% ccfB1B1 = CalcCCF_FFT(SignBarkerB1,SignBarkerB1, 0);
% figure, plot(ccfB1B1);
% title('ccfB1B1 periodic');
% 
% ccfB2B2 = CalcCCF_FFT(SignBarkerB2,SignBarkerB2, 0);
% figure, plot(ccfB2B2);
% title('ccfB2B2 periodic');
% 
% ccfB1B2 = CalcCCF_FFT(SignBarkerB1,SignBarkerB2, 0);
% figure, plot(ccfB1B2);
% title('ccfB1B2 periodic');


period = 2048;
nInfBits = 2048*4+11;
signalInf_b = 2*randi([0,1],nInfBits,1)-1; %information signal = noise
% signalInf_b = 1:14;
% SignSyncB2      = [-2 -2];
% SignSyncB1      = [-1 -1 -1];
signalInf_B2    = InsertSyncB2(signalInf_b,3*SignBarkerB2, period);
signalInf_B2B1  = InsertSyncB1(signalInf_B2,2*SignBarkerB1, 1000);

figure, plot(signalInf_B2B1);
title('signalInf B2B1');

ccfInfB1 = CalcCCF_FFT(signalInf_B2B1,SignBarkerB1, 0);
figure, plot(ccfInfB1);
title('ccfInfB1');

ccfInfB2 = CalcCCF_FFT(signalInf_B2B1,SignBarkerB2, 0);
figure, plot(ccfInfB2);
title('ccfInfB2');

x = 1:length(ccfInfB1);
figure, plot(x,ccfInfB1,'r',x,ccfInfB2,'b');
title('ccfInfB1 (red) and ccfInfB2 (blue)');

%**************************************************
%***************Receiver***************************
%**************************************************
SignSyncB1  = ccfInfB1;
SignSyncB2  = ccfInfB2;
EstSignal   = signalInf_B2B1;
nInfBits    = nInfBits;
period      = period;

[MaxSignSyncB1,ImaxSignSyncB1] = max(SignSyncB1);  %largest element index
[MinSignSyncB1,IminSignSyncB1] = min(SignSyncB1);  %smalles element index

iAB1 = ImaxSignSyncB1 + length(SignBarkerB1);
iBB1 = IminSignSyncB1-1;
n_packets = fix(nInfBits/period);
packet = zeros(period,n_packets);
lenghth_packet_tail = nInfBits - n_packets*period;  %last packet size if abs(nInfBits/period - fix(nInfBits/period)) > 0
packet_tail = zeros(lenghth_packet_tail,1); %last packet if abs(nInfBits/period - fix(nInfBits/period)) > 0
EstSignal_b = zeros(n_packets*period + lenghth_packet_tail,1);

MaxSignSyncB2     = zeros(n_packets,1);   %MaxSignSync is max(CCFB2)
%MinSignSyncB2     = zeros(n_packets,1);   %MinSignSync is min(CCFB2)
ImaxSignSyncB2     = zeros(n_packets,1);   %iMaxSignSync is imax(CCFB2)
%IminSignSyncB2     = zeros(n_packets,1);   %iMinSignSync is imin(CCFB2)



for i=1:n_packets
    mask = zeros(length(SignSyncB2),1);
    iB = iAB1+i*(period+length(SignBarkerB2));
    iA = iB - 2*length(SignBarkerB2);
    mask(iA:iB) = 1;
    
    [MaxSignSyncB2(i),ImaxSignSyncB2(i)] = max(SignSyncB2.*mask);  %largest element index
    
    x = 1:length(SignSyncB2);
    figure, plot(x,SignSyncB2,x,mask);
end

SignSyncB1_test = zeros(length(SignSyncB1),1);
SignSyncB1_test(iAB1:iBB1) = 1;
x = 1:length(SignSyncB1);
figure, plot(x,SignSyncB1,x,SignSyncB1_test);

if n_packets >0
    packet(:,1) = EstSignal(iAB1:ImaxSignSyncB2(1)-1);
    for i=2:n_packets
        iB = ImaxSignSyncB2(i)-1;
        iA = ImaxSignSyncB2(i-1) + length(SignBarkerB2);
        packet(:,i) = EstSignal(iA:iB);
    end
    ind = 1;
    for i=1:n_packets
        EstSignal_b(ind:ind+period-1) = packet(:,i);
        ind = ind + period;
    end
end

% if (lenghth_packet_tail >0) && (n_packets >0)
%     packet_tail(:) = EstSignal(ImaxSignSyncB2(n_packets) + length(SignBarkerB2):iBB1);
%     EstSignal_b(length(EstSignal_b) - lenghth_packet_tail+1:length(EstSignal_b)) = packet_tail;
% end
% 
% if (lenghth_packet_tail >0) && (n_packets == 0)
%     packet_tail(:) = EstSignal(iAB1:iBB1);
%     EstSignal_b(length(EstSignal_b) - lenghth_packet_tail+1:length(EstSignal_b)) = packet_tail;
% end

if lenghth_packet_tail >0
    if n_packets >0
        packet_tail(:) = EstSignal(ImaxSignSyncB2(n_packets) + length(SignBarkerB2):iBB1);
    else
        packet_tail(:) = EstSignal(iAB1:iBB1);
    end
    EstSignal_b(length(EstSignal_b) - lenghth_packet_tail+1:length(EstSignal_b)) = packet_tail;
end


% SignSyncB2_test = zeros(length(SignSyncB2),1);
% SignSyncB2_test(iAB1+period-length(SignBarkerB2):iAB1+period+length(SignBarkerB2)) = 1;
% x = 1:length(SignSyncB2);
% figure, plot(x,SignSyncB2,x,SignSyncB2_test);
% 
% SignSyncB2_test = zeros(length(SignSyncB2),1);
% SignSyncB2_test(iAB1+2*period:iAB1+2*(period+length(SignBarkerB2))) = 1;
% x = 1:length(SignSyncB2);
% figure, plot(x,SignSyncB2,x,SignSyncB2_test);
% 
% SignSyncB2_test = zeros(length(SignSyncB2),1);
% SignSyncB2_test(iAB1+3*period+length(SignBarkerB2):iAB1+3*(period+length(SignBarkerB2))) = 1;
% x = 1:length(SignSyncB2);
% figure, plot(x,SignSyncB2,x,SignSyncB2_test);
% 
% SignSyncB2_test = zeros(length(SignSyncB2),1);
% SignSyncB2_test(iAB1+4*period+2*length(SignBarkerB2):iAB1+4*(period+length(SignBarkerB2))) = 1;
% x = 1:length(SignSyncB2);
% figure, plot(x,SignSyncB2,x,SignSyncB2_test);
% 
% 
% SignSyncB2_test = zeros(length(SignSyncB2),1);
% SignSyncB2_test(iAB1+period:iBB1) = 1;
% x = 1:length(SignSyncB2);
% figure, plot(x,SignSyncB2,x,SignSyncB2_test);

% Samples = 1;
% nMax = round((iB-iA)/Samples);
% delta = iB-iA;




% iA = iA + fix(Samples/2);
% %EstSignalDecimation = zeros(length(EstSignal),1);
% EstSignal_b =  zeros(nMax,1);
% i = 1;
% for n = iA:Samples:iB
% %    EstSignalDecimation(n) = EstSignal(n);
%     EstSignal_b(i) = EstSignal(n);
%     i = i + 1;
% end
% %****syncronization stop*******
