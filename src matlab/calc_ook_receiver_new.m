% 2019-11-17
% noncoherent reception, information signal estimation, BER calculation
function [est_signal_b, index_a, index_b] = calc_ook_receiver_new(z, samples, F, Fs, sign_barker_b1_long, sign_barker_b2_long, n_inf_bits, period, signal_inf_bits)

%****noncoherent reception start *******
signal_complex = CalcNoncoherentReceptionNew(z, samples, F, Fs);            % signal_complex - complex signal
corr_integral = real(signal_complex).^2 + imag(signal_complex).^2;          % detected amplitude (amplitude envelope quadrature)
%****noncoherent reception stop*******

%****information signal estimation (start)*******
%threshold = 0:0.01:0.6;  %resolver threshold. max(corr_integral)
threshold = linspace(min(corr_integral), max(corr_integral));
n = length(threshold);
BER = -2 * ones(n, 1);
max_sign_sync = -2 * ones(n, 1);
min_sign_sync = -2 * ones(n, 1);

for i = 1:n
    [est_signal_b, max_sign_sync(i), min_sign_sync(i)] = calc_signal_estimation_ook(corr_integral, threshold(i), sign_barker_b1_long, sign_barker_b2_long, samples, n_inf_bits, period, signal_complex); %This function estimates information bits (information signal)
%     if length(est_signal_b) == n_inf_bits
%     %if length(est_signal_b) == length(signal_inf_bits)
%         BER(i) = mean(abs(est_signal_b - signal_inf_bits) / 2);   %The bit error rate (BER)
%     else
%         disp(['Error. Can not calculate BER, because of different array size. length(est_signal_b) = ', num2str(length(est_signal_b)), ',  length(signal_inf_bits) = ', num2str(length(signal_inf_bits))]);
%     end
    BER(i) = calc_ber(signal_inf_bits, est_signal_b, n_inf_bits);
end
%****information signal estimation (stop)*******

m = -2;
i = -2;
[m i] = max(max_sign_sync - min_sign_sync);
disp(['threshold = ', num2str(threshold(i))]);
disp(['max_sign_sync-min_sign_sync = ', num2str(m)]);
disp(['BER = ', num2str(BER(i))]);
[signal_complex] = CalcNoncoherentReceptionNew(z, samples, F, Fs);      %signal_complex - complex signal
corr_integral = real(signal_complex).^2 + imag(signal_complex).^2;       %detected amplitude (amplitude envelope quadrature)
[est_signal_b, ~, ~, ~, signal_constel, index_a, index_b] = calc_signal_estimation_ook(corr_integral, threshold(i), sign_barker_b1_long, sign_barker_b2_long, samples, n_inf_bits, period, signal_complex); % this function estimates information bits (information signal)
thr = threshold(i);

plot_time(corr_integral, Fs, 'sec', 'corr integral')

x = 1:length(z);
x = x / Fs;
figure, plot(x, corr_integral, 'r', x, z, 'b');
xlabel('sec');
title('SignAmp (r) and z (b)');

c = linspace(1, 10, length(signal_constel));                         % from black to yellow
figure, scatter(real(signal_constel), imag(signal_constel), [], c);     % create a scatter plot and vary the circle color.
hold on;
theta = linspace(0, 2 * pi);
r = sqrt(thr);
x = r * cos(theta);
y = r * sin(theta);
plot(x, y);
axis equal; %Use the same length for the data units along each axis.
xlabel('In Phase');
ylabel('Quadrature');
title('Signal Constellation');