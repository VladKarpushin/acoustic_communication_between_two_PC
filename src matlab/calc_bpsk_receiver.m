% 2019-12-15
% coherent reception, information signal estimation, BER calculation

function [est_signal_b, index_a, index_b] = calc_bpsk_receiver(z, samples, F, Fs, sign_barker_long, n_inf_bits, signal_inf_bits)

%*******PLL start ******
n = 2 * ceil(samples / 4);          %number of PLL iterations n = Pi/4
PLL_offset_n = 0:n - 1;           %PLL_offset_n - %sin offset (used for PLL). It is equal to PLL_offset_n  = ceil(samples/4). PLL_offset_n = 0 means there is not offset

max_abs_corr_integral   = -2 * ones(n, 1);   %max(CCF received signal and sin wave), max(correlation integral)
BER                     = -2 * ones(n, 1);   %BER is bit error rate
max_sign_sync           = -2 * ones(n, 1);   %max_sign_sync is max(CCF)
min_sign_sync           = -2 * ones(n, 1);   %min_sign_sync is min(CCF)
delta                   = -2 * ones(n, 1);   %delta is difference between index(min_sign_sync) and index(max_sign_sync)
std_sign_sync           = -2 * ones(n, 1);   %std_sign_sync is std(CCF)
threshold = 0;                      %resolver threshold. Should be zero for BPSK

for i = 1:n
    [corr_integral, tmp] = CalcCoherentReceptionNew3(z, samples, F, Fs, PLL_offset_n(i));   %coherent reception
    [est_signal_b, max_sign_sync(i), min_sign_sync(i), Err delta(i), std_sign_sync(i)] = CalcSignalEstimationNew4(corr_integral, threshold, sign_barker_long, samples, tmp); %This function estimates information bits (information signal)
    max_abs_corr_integral(i) = max(abs(corr_integral));
    %if length(est_signal_b) == length(signal_inf_bits)                  %check size
%     if length(est_signal_b) == n_inf_bits
%         BER(i) = mean(abs(est_signal_b - signal_inf_bits) / 2);   %The bit error rate (BER) calculation
%     else
%         disp(['Error. Can not calculate BER, because of different array size. length(est_signal_b) = ', num2str(length(est_signal_b)), ',  length(signal_inf_bits) = ', num2str(length(signal_inf_bits))]);
%     end
    BER(i) = calc_ber(signal_inf_bits, est_signal_b, n_inf_bits);

end
ErrSyst = n_inf_bits * samples - delta; %systematic error between n_inf_bits*samples and delta
%PLL_offset_vs_BER = [PLL_offset_n' max_abs_corr_integral max_sign_sync min_sign_sync max_sign_sync-min_sign_sync BER delta ErrSyst std_sign_sync];
%*******PLL stop ******

%*******output result (start)*********
m = -2;
i = -2;
[m i] = max(max_sign_sync - min_sign_sync);   %PLL criterion1: max(CCF) - min(CCF) = max
if abs(ErrSyst(i)) > samples / 2      
    [m i] = min(abs(ErrSyst));          %PLL criterion2: systematic error = min
end

disp(['PLL_offset_n = ', num2str(PLL_offset_n(i))]);
disp(['BER = ', num2str(BER(i))]);
disp(['max_sign_sync - min_sign_sync = ', num2str(max_sign_sync(i) - min_sign_sync(i))]);
%disp(['delta = ', num2str(delta(i))]);
%disp(['ErrSyst (systematic error) = ', num2str(ErrSyst(i))]);
%disp(['std_sign_sync = ', num2str(std_sign_sync(i))]);

[corr_integral signal_complex] = CalcCoherentReceptionNew3(z, samples, F, Fs, PLL_offset_n(i));   %coherent reception
[est_signal_b max_sign_sync min_sign_sync Err delta std_sign_sync signal_constel index_a index_b] = CalcSignalEstimationNew4(corr_integral, threshold, sign_barker_long, samples, signal_complex); %This function estimates information bits (information signal)

plot_time(corr_integral, Fs, 'sec', 'corr integral')

x = 1:length(z);
x = x / Fs;
figure, plot(x, corr_integral, 'r', x, z, 'b');
xlabel('sec');
title('corr-integral (r) and z (b)');

c = linspace(1, 10, length(signal_constel));                         % from black to yellow
figure, scatter(real(signal_constel), imag(signal_constel), [], c);     % create a scatter plot and vary the circle color.
axis equal; %Use the same length for the data units along each axis.
xlabel('In Phase');
ylabel('Quadrature');
title('Signal Constellation');