% 2019-12-08
% to display signal parameters

function show_sign_para(kt, F, Fs, n_total_bits, n_inf_bits)

t = kt * n_total_bits / F;   %common transmit time
disp(['Sampling rate = ', num2str(Fs),' Hz']);
disp(['Freq of signal1 = ', num2str(F),' Hz']);
disp(['Duration of one symbol = ', num2str(1000 * kt / F),' ms']);
disp(['BW = ', num2str(2 / (kt / F)),' Hz']);
disp(['Number of total transmitted bits = ', num2str(n_total_bits),' [bits]']);
disp(['Number of information transmitted bits = ', num2str(n_inf_bits),' [bits]']);
disp(['Physical layer gross bitrate = ', num2str(F / kt),' [bits/s]']);
disp(['Net bitrate = ', num2str(fix(n_inf_bits / t)),' [bits/s]']);
disp(['Symbol rate = ', num2str(F / kt),' [Hz]']);
disp(['common transmit time = ', num2str(t),' [s]']);