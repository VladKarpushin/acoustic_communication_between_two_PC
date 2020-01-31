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

est_signal_long = zeros(length(corr_integral), 1);
est_signal_long = (2 * (corr_integral > threshold)) - 1;    %resolver

%**************************************************************
sync_b1 = CalcCCF_FFT(est_signal_long, sign_barker_b1_long, 0);

[max_sync_b1, ind_max_sync_b1] = max(sync_b1);  %largest element index
[min_sync_b1, ind_min_sync_b1] = min(sync_b1);  %smalles element index

ind_a = ind_max_sync_b1 + length(sign_barker_b1_long);
ind_b = ind_min_sync_b1 - 1;

if ind_a > ind_b                                      %inversion check
    ind_a = ind_min_sync_b1 + length(sign_barker_b1_long);
    ind_b = ind_max_sync_b1 - 1;
    est_signal_long = - est_signal_long;
end
%StdSignSync = std(sync_b1(ind_a:ind_b-length(sign_barker_b1_long)));  %std(CCF) between two mainlobes
%delta = ind_b - ind_a;

sync_b2 = CalcCCF_FFT(est_signal_long, sign_barker_b2_long, 0);
% figure, plot(sync_b2);
% title('sync_b2');
% 
% x = 1:length(sync_b1);
% figure, plot(x,sync_b1,'r',x,sync_b2,'b');
% title('sync_b1 (red) and sync_b2 (blue)');

n_packets = fix(n_inf_bits/period);
packet = zeros(period,n_packets);
packet_SignalContell = zeros(period,n_packets);
lenghth_packet_tail = n_inf_bits - n_packets*period;  %last packet size if abs(n_inf_bits/period - fix(n_inf_bits/period)) > 0
packet_tail = zeros(lenghth_packet_tail,1); %last packet if abs(n_inf_bits/period - fix(n_inf_bits/period)) > 0
packet_tail_SignalContell = zeros(lenghth_packet_tail,1); %last packet if abs(n_inf_bits/period - fix(n_inf_bits/period)) > 0
est_signal_b = zeros(n_packets*period + lenghth_packet_tail,1);
signal_constel = zeros(n_packets*period + lenghth_packet_tail,1);

max_sync_b2     = zeros(n_packets,1);   %MaxSignSync is max(CCFB2)
ind_max_sync_b2     = zeros(n_packets,1);   %iMaxSignSync is imax(CCFB2)

for i=1:n_packets
    mask = zeros(length(sync_b2),1);
    iB = ind_a+i*(period*samples+length(sign_barker_b2_long));
    iA = iB - 2*length(sign_barker_b2_long);
    
    if iB <= length(mask)
        mask(iA:iB) = 1;
        tmp = sync_b2.*mask;
        [max_sync_b2(i), ind_max_sync_b2(i)] = max(sync_b2 .* mask);  %largest element index
    else
        %disp(['Error. Wrong iB value) = ',num2str(iB)]);
        Err = 1;
        return;
    end
    
    
%     x = 1:length(sync_b2);
%     figure, plot(x,sync_b2,x,mask);
end


if n_packets >0
    packetLong = est_signal_long(ind_a:ind_max_sync_b2(1)-1);
    packetLong_SignalContell = signal_complex(ind_a:ind_max_sync_b2(1)-1);
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
        iB = ind_max_sync_b2(i)-1;
        iA = ind_max_sync_b2(i-1) + length(sign_barker_b2_long);
        packetLong  = est_signal_long(iA:iB);
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
        packetLong = est_signal_long(ind_max_sync_b2(n_packets) + length(sign_barker_b2_long):ind_b);
        packetLong_SignalContell = signal_complex(ind_max_sync_b2(n_packets) + length(sign_barker_b2_long):ind_b);
    else
        packetLong = est_signal_long(ind_a:ind_b);
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
