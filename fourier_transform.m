%% TRANSFORMADA DE FOURIER - Representaciones 2D en el dominio frecuencial
% Acepta un archivo de audio como señal de entrada.
% Si no se especifica, usa una señal sintética de ejemplo.

clear; clc; close all;

%% ── 0. ARCHIVO DE AUDIO (opcional) ──────────────────────────────────────────
% Indica la ruta al archivo de audio (WAV, MP3, FLAC, OGG, etc.).
% Deja la cadena vacía ('') para usar la señal sintética de ejemplo.
audio_file = '';   % Ejemplo: audio_file = 'mi_audio.wav';

%% ── 1. CARGA O GENERACIÓN DE LA SEÑAL ───────────────────────────────────────
if ~isempty(audio_file) && isfile(audio_file)
    fprintf('Cargando archivo de audio: %s\n', audio_file);
    [signal_raw, fs] = audioread(audio_file);
    % Convertir estéreo (o multicanal) a mono promediando canales
    if size(signal_raw, 2) > 1
        signal_raw = mean(signal_raw, 2);
    end
    signal   = signal_raw(:)';          % vector fila
    duration = length(signal) / fs;
    t        = (0 : length(signal)-1) / fs;
    fprintf('  fs = %d Hz  |  Duración = %.2f s  |  Muestras = %d\n', ...
            fs, duration, length(signal));
elseif ~isempty(audio_file)
    warning('Archivo no encontrado: "%s"\nUsando señal sintética de ejemplo.', audio_file);
    fs = 1000;  duration = 1.0;
    t  = 0 : 1/fs : duration - 1/fs;
    signal = 1.00*sin(2*pi*50*t) + 0.50*sin(2*pi*120*t) + 0.25*sin(2*pi*300*t);
else
    fprintf('No se especificó archivo de audio. Usando señal sintética de ejemplo.\n');
    fs = 1000;  duration = 1.0;
    t  = 0 : 1/fs : duration - 1/fs;
    signal = 1.00*sin(2*pi*50*t) + 0.50*sin(2*pi*120*t) + 0.25*sin(2*pi*300*t);
end

N = length(signal);

%% ── 2. FFT ───────────────────────────────────────────────────────────────────
Y         = fft(signal);
freqs     = (0 : N/2) * fs / N;          % eje de frecuencias (espectro unilateral)
magnitude = abs(Y(1:N/2+1)) / N * 2;     % normalización (x2 por espectro unilateral)
magnitude(1)   = magnitude(1) / 2;       % componente DC sin x2
magnitude(end) = magnitude(end) / 2;     % componente Nyquist sin x2
phase     = angle(Y(1:N/2+1));           % fase (radianes)

%% ── 3. ESPECTROGRAMA (Short-Time Fourier Transform) ─────────────────────────
win_len  = min(256, floor(N / 4));        % ventana adaptada a la longitud de la señal
hop      = floor(win_len / 4);
win      = hann(win_len);
[S, F, T_spec] = spectrogram(signal, win, win_len - hop, win_len, fs);

%% ── 4. FIGURA 1: Señal temporal + Espectro de magnitud + Fase ───────────────
f_max = fs / 2;                            % máxima frecuencia a mostrar

figure('Name','Transformada de Fourier', 'NumberTitle','off', ...
       'Position',[100 100 1100 700]);

% --- 4a. Señal en el tiempo
subplot(3,2,[1 2]);
plot(t, signal, 'b', 'LineWidth', 1.2);
title('Señal en el dominio temporal');
xlabel('Tiempo (s)');  ylabel('Amplitud');
xlim([0 duration]);
grid on;

% --- 4b. Espectro de magnitud (escala lineal)
subplot(3,2,3);
plot(freqs, magnitude, 'r', 'LineWidth', 0.8);
title('Espectro de magnitud (lineal)');
xlabel('Frecuencia (Hz)');  ylabel('|X(f)|');
xlim([0 f_max]);
grid on;

% --- 4c. Espectro de magnitud (escala logarítmica dB)
subplot(3,2,4);
plot(freqs, 20*log10(magnitude + 1e-12), 'm', 'LineWidth', 1.2);
title('Espectro de magnitud (dB)');
xlabel('Frecuencia (Hz)');  ylabel('Magnitud (dB)');
xlim([0 f_max]);  ylim([-80 10]);
grid on;

% --- 4d. Espectro de fase
subplot(3,2,5);
plot(freqs, rad2deg(phase), 'g', 'LineWidth', 1.0);
title('Espectro de fase');
xlabel('Frecuencia (Hz)');  ylabel('Fase (°)');
xlim([0 f_max]);
grid on;

% --- 4e. Espectrograma 2D (tiempo × frecuencia)
subplot(3,2,6);
imagesc(T_spec, F, 20*log10(abs(S) + 1e-12));
axis xy;
colormap(gca, 'jet');
colorbar;
title('Espectrograma 2D (STFT)');
xlabel('Tiempo (s)');  ylabel('Frecuencia (Hz)');
ylim([0 f_max]);
clim([-80 0]);

sgtitle(sprintf('Análisis de Fourier — fs = %d Hz', fs), 'FontSize', 14, 'FontWeight', 'bold');

%% ── 5. FIGURA 2: Representación 2D de la energía espectral ──────────────────
figure('Name','Densidad espectral de potencia 2D', 'NumberTitle','off', ...
       'Position',[150 150 900 500]);

win_pwelch = min(256, floor(N / 4));
[pxx, fpxx] = pwelch(signal, hann(win_pwelch), floor(win_pwelch/2), win_pwelch*2, fs);

subplot(1,2,1);
area(fpxx, 10*log10(pxx), 'FaceColor',[0.2 0.5 0.8], 'FaceAlpha',0.6, 'EdgeColor','b');
title('Densidad Espectral de Potencia (Welch)');
xlabel('Frecuencia (Hz)');  ylabel('PSD (dB/Hz)');
xlim([0 f_max]);
grid on;

% Mapa de calor 2D: ventanas × bins de frecuencia
subplot(1,2,2);
F_mask = F <= f_max;
pcolor(T_spec, F(F_mask), abs(S(F_mask,:)));
shading interp;
colormap(gca, 'hot');
colorbar;
title('Mapa de calor espectral 2D');
xlabel('Tiempo (s)');  ylabel('Frecuencia (Hz)');

sgtitle('Análisis de potencia espectral', 'FontSize', 13, 'FontWeight', 'bold');

%% ── 6. RESUMEN EN CONSOLA ────────────────────────────────────────────────────
fprintf('\n=== Frecuencias dominantes ===\n');
fprintf('%-20s %-12s\n', 'Frecuencia (Hz)', 'Magnitud');
fprintf('%s\n', repmat('-', 1, 34));

threshold = max(magnitude) * 0.05;   % 5% del pico máximo
dom_idx   = find(magnitude > threshold);
[~, ord]  = sort(magnitude(dom_idx), 'descend');
dom_idx   = dom_idx(ord);

for k = 1:min(10, length(dom_idx))
    fprintf('%-20.1f %-12.4f\n', freqs(dom_idx(k)), magnitude(dom_idx(k)));
end
