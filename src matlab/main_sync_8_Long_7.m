%2016-08-06
%sync test with Barker codes of different length
%2016-08-13. added InsertSyncB2 and InsertSyncB1
%2016-08-20 recoverred inf signal

close all,clc,clear all;

%******************************
%*******Transmitter************
%******************************

Samples  =30;
delay = 1000;
nSignBarker = 75;   %quantity of Barker codes in a set.
SignBarkerOneB1 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1]'; %Barker code N=13. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems. Barker codes have length at most 13 and have low correlation sidelobes
SignBarkerOneB2 = [1 1 1 -1 -1 -1 1 -1 -1 1 -1]'; %Barker code N=11. Barker codes, which are subsets of PN sequences, are commonly used for frame synchronization in digital communication systems.

SignBarkerB1 = GetPeriodicBarkerCode(SignBarkerOneB1, nSignBarker);
SignBarkerB2 = GetPeriodicBarkerCode(SignBarkerOneB2, nSignBarker);


period = 2048;
nInfBits = 2048*4+100;
signalInf_b = 2*randi([0,1],nInfBits,1)-1; %information signal = noise

signalInf_B2    = InsertSyncB2(signalInf_b,SignBarkerB2, period);
signalInf_B2B1  = InsertSyncB1(signalInf_B2,SignBarkerB1, delay);

%**************************************************
%***************Receiver***************************
%**************************************************

nInfBits    = nInfBits;
period      = period;
Samples     = Samples;

EstSignal           = - Short2Long(signalInf_B2B1, Samples);
SignBarkerB1Long    = Short2Long(SignBarkerB1, Samples);
SignBarkerB2Long    = Short2Long(SignBarkerB2, Samples);

figure, plot(EstSignal);
title('EstSignal (signalInf B2B1)');

SignSyncB1 = CalcCCF_FFT(EstSignal,SignBarkerB1Long, 0);
figure, plot(SignSyncB1);
title('SignSyncB1');

[MaxSignSyncB1,ImaxSignSyncB1] = max(SignSyncB1);  %largest element index
[MinSignSyncB1,IminSignSyncB1] = min(SignSyncB1);  %smalles element index

iAB1 = ImaxSignSyncB1 + length(SignBarkerB1Long);
iBB1 = IminSignSyncB1-1;

if iAB1 > iBB1                  
    iAB1 = IminSignSyncB1 + length(SignBarkerB1Long);
    iBB1 = ImaxSignSyncB1-1;
    EstSignal = -EstSignal;
end

SignSyncB2 = CalcCCF_FFT(EstSignal,SignBarkerB2Long, 0);
figure, plot(SignSyncB2);
title('SignSyncB2');

x = 1:length(SignSyncB1);
figure, plot(x,SignSyncB1,'r',x,SignSyncB2,'b');
title('SignSyncB1 (red) and SignSyncB2 (blue)');




n_packets = fix(nInfBits/period);
packet = zeros(period,n_packets);
lenghth_packet_tail = nInfBits - n_packets*period;  %last packet size if abs(nInfBits/period - fix(nInfBits/period)) > 0
packet_tail = zeros(lenghth_packet_tail,1); %last packet if abs(nInfBits/period - fix(nInfBits/period)) > 0
EstSignal_b = zeros(n_packets*period + lenghth_packet_tail,1);

MaxSignSyncB2     = zeros(n_packets,1);   %MaxSignSync is max(CCFB2)
ImaxSignSyncB2     = zeros(n_packets,1);   %iMaxSignSync is imax(CCFB2)

for i=1:n_packets
    mask = zeros(length(SignSyncB2),1);
    iB = iAB1+i*(period*Samples+length(SignBarkerB2Long));
    iA = iB - 2*length(SignBarkerB2Long);
    mask(iA:iB) = 1;
    
    [MaxSignSyncB2(i),ImaxSignSyncB2(i)] = max(SignSyncB2.*mask);  %largest element index
    
    x = 1:length(SignSyncB2);
    figure, plot(x,SignSyncB2,x,mask);
end

% SignSyncB1_test = zeros(length(SignSyncB1),1);
% SignSyncB1_test(iAB1:iBB1) = 1;
% x = 1:length(SignSyncB1);
% figure, plot(x,SignSyncB1,x,SignSyncB1_test);

if n_packets >0
    packetLong = EstSignal(iAB1:ImaxSignSyncB2(1)-1);
    packet(:,1) = Long2Short(packetLong,Samples);
    for i=2:n_packets
        iB = ImaxSignSyncB2(i)-1;
        iA = ImaxSignSyncB2(i-1) + length(SignBarkerB2Long);
        packetLong  = EstSignal(iA:iB);
        packet(:,i) = Long2Short(packetLong,Samples);
    end
    ind = 1;
    for i=1:n_packets
        EstSignal_b(ind:ind+period-1) = packet(:,i);
        ind = ind + period;
    end
end

if lenghth_packet_tail >0
    if n_packets >0
        packetLong = EstSignal(ImaxSignSyncB2(n_packets) + length(SignBarkerB2Long):iBB1);
    else
        packetLong = EstSignal(iAB1:iBB1);
    end
    packet_tail = Long2Short(packetLong,Samples);
    EstSignal_b(length(EstSignal_b) - lenghth_packet_tail+1:length(EstSignal_b)) = packet_tail;
end

if length(EstSignal_b) == length(signalInf_b)                  %check size
    BER = mean(abs(EstSignal_b - signalInf_b)/2);   %The bit error rate (BER) calculation
else
    disp(['Error. Can not calculate BER, because of different array size. length(EstSignal_b) = ',num2str(length(EstSignal_b)), ',  length(signalInf_b) = ', num2str(length(signalInf_b))]);
end
[signalInf_b EstSignal_b]
BER
