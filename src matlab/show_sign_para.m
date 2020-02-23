% 2019-12-08
% to display signal parameters

function show_sign_para(kt, F, Fs, n_total_bits, n_inf_bits, delay, n_total_bits_q)

M  = 2; % number of bits per modulation symbol. 1 - for bpsk and ook, 2 - for qpsk
if nargin < 7
    M = 1;
    n_total_bits_q = 0;
end

disp(['Sampling rate = ', num2str(Fs),' Hz']);
disp(['f0  = carrier frequency = ', num2str(F),' Hz']);
disp(['K = samples per one symbol = ', num2str(kt * Fs / F),' [samples]']);
disp(['BW = ', num2str(2 / (kt / F)),' Hz']);

%disp(['Number of total transmitted bits = ', num2str(M * n_total_bits - delay),' [bits]']); % should be n_total_bits_i + n_total_bits_q - delay
disp(['Number of total transmitted bits = ', num2str(n_total_bits + n_total_bits_q - delay),' [bits]']); % should be n_total_bits_i + n_total_bits_q - delay
disp(['Number of information transmitted bits = ', num2str(n_inf_bits),' [bits]']);

disp(['Duration of one symbol = ', num2str(1000 * kt / F),' ms']);
disp(['Symbol rate = baud rate = ', num2str(F / kt),' [Hz]']);
disp(['Physical layer gross bitrate = ', num2str(M * F / kt),' [bits / s]']); % For QPSK should be 2 * F / kt

t = kt * (n_total_bits - delay) / F;   %common transmit time
disp(['Net bitrate old = ', num2str(fix(n_inf_bits / t)),' [bits / s]']);

disp(['common transmit time = ', num2str(t),' [s]']);