%function [est_signal_b MaxSignSync MinSignSync Err delta StdSignSync] = CalcSignalEstimationNew2B1B2(CorrIntegral,threshold, SignBarkerB1Long,SignBarkerB2Long, Samples,nInfBits,period)
%This function estimates information bits (information signal) from B1B2
%signal
%without plots
%2016-11-27 added code for signal constellation
%2016-12-03 added index_a and index_b export for SNR calculation

function [est_signal_b, max_sync_b1, min_sync_b1, Err, signal_constel, ind_a, ind_b] = calc_signal_estimation_ook(corr_integral, threshold, sign_barker_b1_long, sign_barker_b2_long, samples, n_inf_bits, period, signal_complex)
% input:
% 	corr_integral    - cross correlation function (CCF) received signal SignR and sin wave. Another name is correlation integral
% 	threshold       - resolver threshold. Should be zero for BPSK
%   SignBarkerLong  -    sync signal(long)
%   samples         - quantity of samples per one symbol
%   signal_complex   - complex signal
% output:
% 	est_signal_b     - estimated information bits
%   Err             - error information
%   MaxSignSync     - MaxSignSync is max(CCF)
%   MinSignSync     - MinSignSync is min(CCF)
%   delta           - delta is difference between index(MinSignSync) and index(MaxSignSync)
%   StdSignSync     - StdSignSync is std(CCF) between two mainlobes
%   signal_constel   - signal constellation
%   ind_a            - index of first symbol. For SNR estimation
%   ind_b            - index of last symbol. For SNR estimation

Err = 0;
max_sync_b1 = 0;
min_sync_b1 = 0;
est_signal_b = 0;
signal_constel = 0;
ind_a = 0;
ind_b = 0;

EstSignal = zeros(length(corr_integral),1);
EstSignal = (2*(corr_integral > threshold))-1;    %resolver

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
% est_signal_b = Long2Short(EstSignal(iA:iB),samples);
%delta = iB-iA;
%****syncronization stop*******


%**************************************************************
%new part
SignSyncB1 = CalcCCF_FFT(EstSignal,sign_barker_b1_long, 0);

[max_sync_b1,ImaxSignSyncB1] = max(SignSyncB1);  %largest element index
[min_sync_b1,IminSignSyncB1] = min(SignSyncB1);  %smalles element index

ind_a = ImaxSignSyncB1 + length(sign_barker_b1_long);
ind_b = IminSignSyncB1-1;

if ind_a > ind_b                                      %inversion check
    ind_a = IminSignSyncB1 + length(sign_barker_b1_long);
    ind_b = ImaxSignSyncB1-1;
    EstSignal = -EstSignal;
end
StdSignSync = std(SignSyncB1(ind_a:ind_b-length(sign_barker_b1_long)));  %std(CCF) between two mainlobes
delta = ind_b-ind_a;

SignSyncB2 = CalcCCF_FFT(EstSignal,sign_barker_b2_long, 0);
% figure, plot(SignSyncB2);
% title('SignSyncB2');
% 
% x = 1:length(SignSyncB1);
% figure, plot(x,SignSyncB1,'r',x,SignSyncB2,'b');
% title('SignSyncB1 (red) and SignSyncB2 (blue)');

n_packets = fix(n_inf_bits/period);
packet = zeros(period,n_packets);
packet_SignalContell = zeros(period,n_packets);
lenghth_packet_tail = n_inf_bits - n_packets*period;  %last packet size if abs(n_inf_bits/period - fix(n_inf_bits/period)) > 0
packet_tail = zeros(lenghth_packet_tail,1); %last packet if abs(n_inf_bits/period - fix(n_inf_bits/period)) > 0
packet_tail_SignalContell = zeros(lenghth_packet_tail,1); %last packet if abs(n_inf_bits/period - fix(n_inf_bits/period)) > 0
est_signal_b = zeros(n_packets*period + lenghth_packet_tail,1);
signal_constel = zeros(n_packets*period + lenghth_packet_tail,1);

MaxSignSyncB2     = zeros(n_packets,1);   %MaxSignSync is max(CCFB2)
ImaxSignSyncB2     = zeros(n_packets,1);   %iMaxSignSync is imax(CCFB2)

for i=1:n_packets
    mask = zeros(length(SignSyncB2),1);
    iB = ind_a+i*(period*samples+length(sign_barker_b2_long));
    iA = iB - 2*length(sign_barker_b2_long);
    
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
    packetLong = EstSignal(ind_a:ImaxSignSyncB2(1)-1);
    packetLong_SignalContell = signal_complex(ind_a:ImaxSignSyncB2(1)-1);
    p = Long2Short(packetLong,samples);
    p_SignalContell = Long2Short(packetLong_SignalContell,samples);
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
        iA = ImaxSignSyncB2(i-1) + length(sign_barker_b2_long);
        packetLong  = EstSignal(iA:iB);
        packetLong_SignalContell = signal_complex(iA:iB);
        p = Long2Short(packetLong,samples);
        p_SignalContell = Long2Short(packetLong_SignalContell,samples);
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
        est_signal_b(ind:ind+period-1) = packet(:,i);
        signal_constel(ind:ind+period-1) = packet_SignalContell(:,i);
        ind = ind + period;
    end
end

if lenghth_packet_tail >0
    if n_packets >0
        packetLong = EstSignal(ImaxSignSyncB2(n_packets) + length(sign_barker_b2_long):ind_b);
        packetLong_SignalContell = signal_complex(ImaxSignSyncB2(n_packets) + length(sign_barker_b2_long):ind_b);
    else
        packetLong = EstSignal(ind_a:ind_b);
        packetLong_SignalContell = signal_complex(ind_a:ind_b);
    end
    p = Long2Short(packetLong,samples);
    p_SignalContell = Long2Short(packetLong_SignalContell,samples);
    if length(p) == lenghth_packet_tail %there was mistake here (period), rectified 20161015
        packet_tail = p;
        packet_tail_SignalContell = p_SignalContell;
    else
        %disp(['Error. Packet size is wrong) = ',num2str(length(p))]);
        Err = 1;
        return;
    end
    %packet_tail = Long2Short(packetLong,samples);
    est_signal_b(length(est_signal_b) - lenghth_packet_tail+1:length(est_signal_b)) = packet_tail;
    signal_constel(length(est_signal_b) - lenghth_packet_tail+1:length(est_signal_b)) = packet_tail_SignalContell;
end


if abs(length(est_signal_b)/8 - fix(length(est_signal_b)/8)) > 0                  %checking if length of est_signal_b is multiple with 8
    disp(['Error. abs(length(est_signal_b)/8 - fix(length(est_signal_b)/8)) > 0. length(est_signal_b) = ',num2str(length(est_signal_b))]);
    Err = 1;
    return;
end
end
