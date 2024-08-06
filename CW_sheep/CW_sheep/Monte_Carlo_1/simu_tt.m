%% harmonic

clear;

fs = 100;  % Sampling freq 100Hz
T = 60;  % total measurement time
t = 0:1/fs:T-1/fs;  

%% ppg generation(time domain)
% MHR, FHR freq
mhr_freq = mean([1.5, 1.75]);  % MHR
fhr_freq = mean([3.5, 3.75]);  % FHR

resp_freq = 0.25;  % Respiration frequency (0.25 Hz corresponds to 15 breaths per minute)

% Generate harmonic frequencies for MHR
harmonic_freqs = mhr_freq * (2:4);  % Second, third, and fourth harmonics

% MHR, FHR, respiration, and harmonic signals
mhr_signal = 0.5 * sin(2 * pi * mhr_freq * t);  % MHR
harmonic_signals = 0.25 * sin(2 * pi * harmonic_freqs(1) * t) + ... % Harmonics
                   0.125 * sin(2 * pi * harmonic_freqs(2) * t) + ...
                   0.0625 * sin(2 * pi * harmonic_freqs(3) * t);
fhr_signal = 0.2 * sin(2 * pi * fhr_freq * t);  % FHR
resp_signal = 0.1 * sin(2 * pi * resp_freq * t);  % Respiration

% DC component
dc_component = 13; % Adjusted DC value

% PPG generation
ppg_signal = mhr_signal + harmonic_signals + fhr_signal + resp_signal + dc_component + 0.05 * randn(size(t));
% PPG plot
figure;
plot(t, ppg_signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Synthetic PPG Signal with DC and Harmonic Components');

% Time domain relative change calculation
Imax = max(ppg_signal);
Imin = min(ppg_signal);
relative_change_time_domain = (Imax - Imin) / Imin;
disp(['Time Domain Relative Change: ', num2str(relative_change_time_domain)]);

%% Frequency domain analysis with STFT
window = hamming(256);  % window function
noverlap = 128;  % overlap number
nfft = 512;  % FFT length

% STFT
[S, F, T] = stft(ppg_signal, fs, 'Window', window, 'OverlapLength', noverlap, 'FFTLength', nfft);
P = abs(S).^2; 

% Find indices for DC, MHR, FHR, and harmonics in the frequency domain
dc_index = F == 0;
mhr_index = F >= mhr_freq & F <= mhr_freq*2;
fhr_index = F >= fhr_freq & F <= fhr_freq*2;
harmonic_indices = arrayfun(@(hf) F >= hf & F <= hf*2, harmonic_freqs, 'UniformOutput', false);
harmonic_indices = any(cat(1, harmonic_indices{:}), 1);  % Combine all harmonic indices


S_dB = 10*log10(P);

% Plot only the positive frequencies of the Spectrogram within 1-10 Hz range
positiveFrequenciesIndex = F >= 0 & F <= 5;
S_dB_positive = S_dB(positiveFrequenciesIndex, :);
F_positive = F(positiveFrequenciesIndex);

% Spectrogram plot
figure;
imagesc(T, F_positive, S_dB_positive);
axis xy; % Flips the y-axis so the lowest frequency is at the bottom
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram with Positive Frequencies (1-5 Hz)');
colorbar;
colormap jet;


% FFT
fft_result = fft(ppg_signal);
N = length(ppg_signal);  
frequencies = (0:N-1)*(fs/N);  
power = abs(fft_result).^2/N;  

% Find indices for DC, MHR, FHR, and respiration in the frequency domain
dc_index = frequencies == 0;
mhr_index = (frequencies >= mhr_freq) & (frequencies <= mhr_freq*1.5);  % Consider only the fundamental frequency
fhr_index = (frequencies >= fhr_freq) & (frequencies <= fhr_freq*1.5);
resp_index = (frequencies >= resre'sp_freq) & (frequencies <= resp_freq*1.5);

% Exclude DC, FHR, and respiration from total power calculation
total_power = sum(power) - power(dc_index) - sum(power(fhr_index)) - sum(power(resp_index));

% Calculate relative change in the frequency domain using FFT (excluding DC, FHR, and respiration)
mhr_power = sum(power(mhr_index));
relative_change_freq_domain_fft = mhr_power / total_power;
disp(['Frequency Domain Relative Change using FFT (excluding DC, FHR, and respiration): ', num2str(relative_change_freq_domain_fft)]);


% Plot only the positive frequencies up to 8 Hz
positiveFrequencies = frequencies(1:floor(N/2)+1);
positivePower = power(1:floor(N/2)+1);

figure;
plot(positiveFrequencies, positivePower);
xlabel('Frequency (Hz)');
ylabel('Power');
title('FFT Result: Power Spectrum within 1-8 Hz Range');
xlim([0.5 8]);  % Set x-axis limits to 1-8 Hz



% Time domain relative change calculation for FHR (Method 1: removing only MHR)
fhr_signal_filtered_1 = ppg_signal - mhr_signal;
Imax_fhr_1 = max(fhr_signal_filtered_1);
Imin_fhr_1 = min(fhr_signal_filtered_1);
relative_change_time_domain_fhr_1 = (Imax_fhr_1 - Imin_fhr_1) / Imin_fhr_1;
disp(['FHR Time Domain Relative Change (Method 1: removing only MHR): ', num2str(relative_change_time_domain_fhr_1)]);

% Time domain relative change calculation for FHR (Method 2: removing MHR and respiration)
fhr_signal_filtered_2 = ppg_signal - mhr_signal - resp_signal;
Imax_fhr_2 = max(fhr_signal_filtered_2);
Imin_fhr_2 = min(fhr_signal_filtered_2);
relative_change_time_domain_fhr_2 = (Imax_fhr_2 - Imin_fhr_2) / Imin_fhr_2;
disp(['FHR Time Domain Relative Change (Method 2: removing MHR and respiration): ', num2str(relative_change_time_domain_fhr_2)]);

% Time domain relative change calculation for FHR (Method 3: removing MHR, respiration, and harmonics)
fhr_signal_filtered_3 = ppg_signal - mhr_signal - resp_signal - harmonic_signals;
Imax_fhr_3 = max(fhr_signal_filtered_3);
Imin_fhr_3 = min(fhr_signal_filtered_3);
relative_change_time_domain_fhr_3 = (Imax_fhr_3 - Imin_fhr_3) / Imin_fhr_3;
disp(['FHR Time Domain Relative Change (Method 3: removing MHR, respiration, and harmonics): ', num2str(relative_change_time_domain_fhr_3)]);

% FFT
fft_result = fft(ppg_signal);
N = length(ppg_signal);  
frequencies = (0:N-1)*(fs/N);  
power = abs(fft_result).^2/N;  

% Find indices for DC, MHR, FHR, and respiration in the frequency domain
dc_index = frequencies == 0;
mhr_index = (frequencies >= mhr_freq) & (frequencies <= mhr_freq*1.5);  % Consider only the fundamental frequency
fhr_index = (frequencies >= fhr_freq) & (frequencies <= fhr_freq*1.5);
resp_index = (frequencies >= resp_freq) & (frequencies <= resp_freq*1.5);

% Exclude DC, FHR, and respiration from total power calculation
total_power = sum(power) - power(dc_index) - sum(power(fhr_index)) - sum(power(resp_index));

% Calculate relative change in the frequency domain using FFT (excluding DC, FHR, and respiration) for MHR
mhr_power = sum(power(mhr_index));
relative_change_freq_domain_fft = mhr_power / total_power;
disp(['MHR Frequency Domain Relative Change using FFT (excluding DC, FHR, and respiration): ', num2str(relative_change_freq_domain_fft)]);

% Frequency domain relative change calculation for FHR
fhr_power = sum(power(fhr_index));
relative_change_freq_domain_fft_fhr = fhr_power / total_power;
disp(['FHR Frequency Domain Relative Change using FFT (excluding DC, MHR, and respiration): ', num2str(relative_change_freq_domain_fft_fhr)]);





