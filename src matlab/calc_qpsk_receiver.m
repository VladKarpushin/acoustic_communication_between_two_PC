% 2019-12-15
% coherent reception, information signal estimation, BER calculation

function [est_signal_b, ind_a, ind_b] = calc_qpsk_receiver(z, samples, F, Fs, sign_barker_b1_long, sign_barker_b2_long, n_inf_bits, signal_inf_bits, criterion)

%*******PLL start ******
n = 50;
%phi = pi * (0:n) / n;
phi = linspace(0, pi, n);
max_abs_corr_integral   = -2 * ones(n, 1);   % max(CCF received signal and sin wave), max(correlation integral)
BER                     = -2 * ones(n, 1);   % BER is bit error rate
max_sync_b1             = -2 * ones(n, 1);   % max(CCF of b1)
max_sync_b2             = -2 * ones(n, 1);   % max(CCF of b2)
delta                   = -2 * ones(n, 1);   %delta is difference between index(min_sync_b2) and index(max_sync_b1)
threshold = 0;                      %resolver threshold. Should be zero for BPSK

%ind_a                   = -2 * ones(n, 1);   %delta is difference between index(min_sync_b2) and index(max_sync_b1)

for i = 1:n
    [signal_complex] = calc_signal_complex(z, samples, F, Fs, phi(i));   %coherent reception
    [est_signal_b, max_sync_b1(i), max_sync_b2(i), ~, ~, ind_a, ind_b] = calc_signal_estimation_bpsk(threshold, sign_barker_b1_long, sign_barker_b2_long, samples, signal_complex, 'qpsk'); %This function estimates information bits (information signal)
    delta(i) = ind_b - ind_a;
    max_abs_corr_integral(i) = max(abs(real(signal_complex)));
    %BER(i) = calc_ber(signal_inf_bits, est_signal_b, n_inf_bits, 'debug');
    BER(i) = calc_ber(signal_inf_bits, est_signal_b, n_inf_bits);
end
err_syst = (n_inf_bits / 2) * samples - delta; %systematic error between n_inf_bits*samples and delta
%PLL_offset_vs_BER = [PLL_offset_n' max_abs_corr_integral max_sync_b1 min_sync_b2 max_sync_b1-min_sync_b2 BER delta err_syst std_sign_sync];
%*******PLL stop ******

%*******output result (start)*********
m = -2;
i = -2;
disp(['phi precision = ', num2str(180 / n), ' [degrees]']);

%if criterion == 'c1'
if strcmp(criterion, 'c1')
    [m i] = max(max_sync_b1 + max_sync_b2);   %PLL criterion1: max(CCF_b1) + max(CCF_b2) = max
    disp(['PLL criterion1: max(max_sync_b1 + max_sync_b2)']);
    disp(['phi c1 = ', num2str(phi(i)), ' [radians]']);
    disp(['phi c1 = ', num2str(phi(i) * 180 / pi), ' [degrees]']);
    disp(['BER c1 = ', num2str(BER(i))]);
%elseif criterion == 'c2'
elseif strcmp(criterion, 'c2')
    [m i] = max(max_sync_b1);   %PLL criterion2: max(max_sync_b1)
    disp(['PLL criterion2: max(max_sync_b1)']);
    disp(['phi c2 = ', num2str(phi(i)), ' [radians]']);
    disp(['phi c2 = ', num2str(phi(i) * 180 / pi), ' [degrees]']);
    disp(['BER c2 = ', num2str(BER(i))]);
else
	disp(['Calc_bpsk_receiver. Wrong criterion']);
end

if abs(err_syst(i)) > samples / 2      
    [m i] = min(abs(err_syst));          %PLL criterion3: systematic error = min
    disp(['PLL criterion3: systematic error = min']);
    disp(['phi c3 = ', num2str(phi(i)), ' [radians]']);
    disp(['phi c3 = ', num2str(phi(i) * 180 / pi), ' [degrees]']);
    disp(['BER c3 = ', num2str(BER(i))]);
end

[signal_complex] = calc_signal_complex(z, samples, F, Fs, phi(i));   %coherent reception
[est_signal_b, ~, ~, ~, signal_constel, ind_a, ind_b] = calc_signal_estimation_bpsk(threshold, sign_barker_b1_long, sign_barker_b2_long, samples, signal_complex, 'qpsk'); %This function estimates information bits (information signal)

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