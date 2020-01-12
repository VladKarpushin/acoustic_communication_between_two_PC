% 2019-12-15
% coherent reception, information signal estimation, BER calculation

function [est_signal_b, ind_a, ind_b] = calc_bpsk_receiver(z, samples, F, Fs, sign_barker_b1_long, sign_barker_b2_long, n_inf_bits, signal_inf_bits)

%*******PLL start ******
n = 20;

phi = pi * (0:n) / n;
max_abs_corr_integral   = -2 * ones(n, 1);   %max(CCF received signal and sin wave), max(correlation integral)
BER                     = -2 * ones(n, 1);   %BER is bit error rate
max_sync_b1             = -2 * ones(n, 1);   %max_sync_b1 is max(CCF)
min_sync_b2             = -2 * ones(n, 1);   %min_sync_b2 is min(CCF)
delta                   = -2 * ones(n, 1);   %delta is difference between index(min_sync_b2) and index(max_sync_b1)
std_sign_sync           = -2 * ones(n, 1);   %std_sign_sync is std(CCF)
threshold = 0;                      %resolver threshold. Should be zero for BPSK

for i = 1:n
    [signal_complex] = calc_coherent_reception_new4(z, samples, F, Fs, phi(i));   %coherent reception
    [est_signal_b, max_sync_b1(i), min_sync_b2(i), Err delta(i), std_sign_sync(i)] = CalcSignalEstimationNew4(threshold, sign_barker_b1_long, sign_barker_b2_long, samples, signal_complex); %This function estimates information bits (information signal)
    
    max_abs_corr_integral(i) = max(abs(real(signal_complex)));
    BER(i) = calc_ber(signal_inf_bits, est_signal_b, n_inf_bits);
end
ErrSyst = n_inf_bits * samples - delta; %systematic error between n_inf_bits*samples and delta
%PLL_offset_vs_BER = [PLL_offset_n' max_abs_corr_integral max_sync_b1 min_sync_b2 max_sync_b1-min_sync_b2 BER delta ErrSyst std_sign_sync];
%*******PLL stop ******

%*******output result (start)*********
m = -2;
i = -2;
[m i] = max(max_sync_b1 + min_sync_b2);   %PLL criterion1: max(CCF) - min(CCF) = max
if abs(ErrSyst(i)) > samples / 2      
    [m i] = min(abs(ErrSyst));          %PLL criterion2: systematic error = min
    disp(['PLL criterion2: systematic error = min']);
end

disp(['phi = ', num2str(phi(i)), ' [radians]']);
disp(['phi = ', num2str(phi(i) * 180 / pi), ' [degrees]']);
disp(['phi precision = ', num2str(180 / n), ' [degrees]']);
disp(['BER = ', num2str(BER(i))]);
disp(['max_sync_b1 - min_sync_b2 = ', num2str(max_sync_b1(i) - min_sync_b2(i))]);

[signal_complex] = calc_coherent_reception_new4(z, samples, F, Fs, phi(i));   %coherent reception
[est_signal_b max_sync_b1 min_sync_b2 Err delta std_sign_sync signal_constel ind_a ind_b] = CalcSignalEstimationNew4(threshold, sign_barker_long, samples, signal_complex); %This function estimates information bits (information signal)

plot_time(real(signal_complex), Fs, 'sec', 'corr integral')

x = 1:length(z);
x = x / Fs;
figure, plot(x, real(signal_complex), 'r', x, z, 'b');
xlabel('sec');
title('corr-integral (r) and z (b)');

c = linspace(1, 10, length(signal_constel));                         % from black to yellow
figure, scatter(real(signal_constel), imag(signal_constel), [], c);     % create a scatter plot and vary the circle color.
axis equal; %Use the same length for the data units along each axis.
xlabel('In Phase');
ylabel('Quadrature');
title('Signal Constellation');