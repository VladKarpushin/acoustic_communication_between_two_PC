%function [EstSignal_b MaxSignSync MinSignSync Err delta StdSignSync] = CalcSignalEstimationNew2B1B2(CorrIntegral,threshold, SignBarkerB1Long,SignBarkerB2Long, Samples,nInfBits,period)
%This function estimates information bits (information signal) from B1B2
%signal
%without plots
%2016-11-27 added code for signal constellation
%2016-12-03 added iAB1 and iBB1 export for SNR calculation

function [EstSignal_b MaxSignSyncB1 MinSignSyncB1 Err delta StdSignSync SignalContell iAB1 iBB1] = CalcSignalEstimationNew4B1B2(CorrIntegral,threshold, SignBarkerB1Long,SignBarkerB2Long, Samples,nInfBits,period,SignalComplex)
% input:
% 	CorrIntegral    - cross correlation function (CCF) received signal SignR and sin wave. Another name is correlation integral
% 	threshold       - resolver threshold. Should be zero for BPSK
%   SignBarkerLong  -    sync signal(long)
%   Samples         - quantity of samples per one symbol
%   SignalComplex   - complex signal
% output:
% 	EstSignal_b     - estimated information bits
%   Err             - error information
%   MaxSignSync     - MaxSignSync is max(CCF)
%   MinSignSync     - MinSignSync is min(CCF)
%   delta           - delta is difference between index(MinSignSync) and index(MaxSignSync)
%   StdSignSync     - StdSignSync is std(CCF) between two mainlobes
%   SignalContell   - signal constellation
%   iAB1            - index of first symbol. For SNR estimation
%   iBB1            - index of last symbol. For SNR estimation

Err = 0;
MaxSignSyncB1 = 0;
MinSignSyncB1 = 0;
EstSignal_b = 0;
SignalContell = 0;
iAB1 = 0;
iBB1 = 0;

EstSignal = zeros(length(CorrIntegral),1);
EstSignal = (2*(CorrIntegral > threshold))-1;    %resolver

%****syncronization start*******
% [SignSync Err] = CalcCCF_FFT(EstSignal,SignBarkerLong,0);
% 
% [MaxSignSync,ImaxSignSync] = max(SignSync);  %largest element index
% [MinSignSync,IminSignSync] = min(SignSync);  %smalles element index

% iA = ImaxSignSync + length(SignBarkerLong);
% iB = IminSignSync-1;
% if iA > iB                  
%     iA = IminSignSync + length(SignBarkerLong);
%     iB = ImaxSignSync-1;
%     EstSignal = -EstSignal;
% end
% StdSignSync = std(SignSync(iA:iB-length(SignBarkerLong)));  %std(CCF) between two mainlobes
% 
% EstSignal_b = Long2Short(EstSignal(iA:iB),Samples);
%delta = iB-iA;
%****syncronization stop*******


%**************************************************************
%new part
SignSyncB1 = CalcCCF_FFT(EstSignal,SignBarkerB1Long, 0);

[MaxSignSyncB1,ImaxSignSyncB1] = max(SignSyncB1);  %largest element index
[MinSignSyncB1,IminSignSyncB1] = min(SignSyncB1);  %smalles element index

iAB1 = ImaxSignSyncB1 + length(SignBarkerB1Long);
iBB1 = IminSignSyncB1-1;

if iAB1 > iBB1                                      %inversion check
    iAB1 = IminSignSyncB1 + length(SignBarkerB1Long);
    iBB1 = ImaxSignSyncB1-1;
    EstSignal = -EstSignal;
end
StdSignSync = std(SignSyncB1(iAB1:iBB1-length(SignBarkerB1Long)));  %std(CCF) between two mainlobes
delta = iBB1-iAB1;

SignSyncB2 = CalcCCF_FFT(EstSignal,SignBarkerB2Long, 0);
% figure, plot(SignSyncB2);
% title('SignSyncB2');
% 
% x = 1:length(SignSyncB1);
% figure, plot(x,SignSyncB1,'r',x,SignSyncB2,'b');
% title('SignSyncB1 (red) and SignSyncB2 (blue)');

n_packets = fix(nInfBits/period);
packet = zeros(period,n_packets);
packet_SignalContell = zeros(period,n_packets);
lenghth_packet_tail = nInfBits - n_packets*period;  %last packet size if abs(nInfBits/period - fix(nInfBits/period)) > 0
packet_tail = zeros(lenghth_packet_tail,1); %last packet if abs(nInfBits/period - fix(nInfBits/period)) > 0
packet_tail_SignalContell = zeros(lenghth_packet_tail,1); %last packet if abs(nInfBits/period - fix(nInfBits/period)) > 0
EstSignal_b = zeros(n_packets*period + lenghth_packet_tail,1);
SignalContell = zeros(n_packets*period + lenghth_packet_tail,1);

MaxSignSyncB2     = zeros(n_packets,1);   %MaxSignSync is max(CCFB2)
ImaxSignSyncB2     = zeros(n_packets,1);   %iMaxSignSync is imax(CCFB2)

for i=1:n_packets
    mask = zeros(length(SignSyncB2),1);
    iB = iAB1+i*(period*Samples+length(SignBarkerB2Long));
    iA = iB - 2*length(SignBarkerB2Long);
    
    if iB <= length(mask)
        mask(iA:iB) = 1;
        tmp = SignSyncB2.*mask;
        [MaxSignSyncB2(i),ImaxSignSyncB2(i)] = max(SignSyncB2.*mask);  %largest element index
    else
        %disp(['Error. Wrong iB value) = ',num2str(iB)]);
        Err = 1;
        return;
    end
    
    
%     x = 1:length(SignSyncB2);
%     figure, plot(x,SignSyncB2,x,mask);
end


if n_packets >0
    packetLong = EstSignal(iAB1:ImaxSignSyncB2(1)-1);
    packetLong_SignalContell = SignalComplex(iAB1:ImaxSignSyncB2(1)-1);
    p = Long2Short(packetLong,Samples);
    p_SignalContell = Long2Short(packetLong_SignalContell,Samples);
    if length(p) == period
        packet(:,1) = p;
        packet_SignalContell(:,1) = p_SignalContell;
    else
        %disp(['Error. Packet size is wrong) = ',num2str(length(p))]);
        Err = 1;
        return;
    end
        
    
    for i=2:n_packets
        iB = ImaxSignSyncB2(i)-1;
        iA = ImaxSignSyncB2(i-1) + length(SignBarkerB2Long);
        packetLong  = EstSignal(iA:iB);
        packetLong_SignalContell = SignalComplex(iA:iB);
        p = Long2Short(packetLong,Samples);
        p_SignalContell = Long2Short(packetLong_SignalContell,Samples);
        if length(p) == period
            packet(:,i) = p;
            packet_SignalContell(:,i) = p_SignalContell;
        else
            %disp(['Error. Packet size is wrong) = ',num2str(length(p))]);
            Err = 1;
            return;
        end
    end
    ind = 1;
    for i=1:n_packets
        EstSignal_b(ind:ind+period-1) = packet(:,i);
        SignalContell(ind:ind+period-1) = packet_SignalContell(:,i);
        ind = ind + period;
    end
end

if lenghth_packet_tail >0
    if n_packets >0
        packetLong = EstSignal(ImaxSignSyncB2(n_packets) + length(SignBarkerB2Long):iBB1);
        packetLong_SignalContell = SignalComplex(ImaxSignSyncB2(n_packets) + length(SignBarkerB2Long):iBB1);
    else
        packetLong = EstSignal(iAB1:iBB1);
        packetLong_SignalContell = SignalComplex(iAB1:iBB1);
    end
    p = Long2Short(packetLong,Samples);
    p_SignalContell = Long2Short(packetLong_SignalContell,Samples);
    if length(p) == lenghth_packet_tail %there was mistake here (period), rectified 20161015
        packet_tail = p;
        packet_tail_SignalContell = p_SignalContell;
    else
        %disp(['Error. Packet size is wrong) = ',num2str(length(p))]);
        Err = 1;
        return;
    end
    %packet_tail = Long2Short(packetLong,Samples);
    EstSignal_b(length(EstSignal_b) - lenghth_packet_tail+1:length(EstSignal_b)) = packet_tail;
    SignalContell(length(EstSignal_b) - lenghth_packet_tail+1:length(EstSignal_b)) = packet_tail_SignalContell;
end


if abs(length(EstSignal_b)/8 - fix(length(EstSignal_b)/8)) > 0                  %checking if length of EstSignal_b is multiple with 8
    disp(['Error. abs(length(EstSignal_b)/8 - fix(length(EstSignal_b)/8)) > 0. length(EstSignal_b) = ',num2str(length(EstSignal_b))]);
    Err = 1;
    return;
end
end
