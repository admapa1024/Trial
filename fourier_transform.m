%% TRANSFORMADA DE FOURIER - Representaciones 2D en el dominio frecuencial
% Señal compuesta por 3 sinusoides con amplitudes y frecuencias distintas

clear; clc; close all;

%% ── 1. PARÁMETROS Y GENERACIÓN DE LA SEÑAL ──────────────────────────────────
fs       = 1000;          % Frecuencia de muestreo (Hz)
duration = 1.0;           % Duración (s)
t        = 0 : 1/fs : duration - 1/fs;
N        = length(t);

% Señal = suma de sinusoides (mismas componentes que el script Python)
f1 = 50;   A1 = 1.00;
f2 = 120;  A2 = 0.50;
f3 = 300;  A3 = 0.25;

signal = A1*sin(2*pi*f1*t) + ...
         A2*sin(2*pi*f2*t) + ...
         A3*sin(2*pi*f3*t);

%% ── 2. FFT ───────────────────────────────────────────────────────────────────
Y         = fft(signal);
freqs     = (0 : N/2) * fs / N;          % eje de frecuencias (espectro unilateral)
magnitude = abs(Y(1:N/2+1)) / N * 2;     % normalización (x2 por espectro unilateral)
magnitude(1)   = magnitude(1) / 2;       % componente DC sin x2
magnitude(end) = magnitude(end) / 2;     % componente Nyquist sin x2
phase     = angle(Y(1:N/2+1));           % fase (radianes)

%% ── 3. ESPECTROGRAMA (Short-Time Fourier Transform) ─────────────────────────
win_len  = 128;           % longitud de ventana (muestras)
hop      = 32;            % salto entre ventanas
win      = hann(win_len); % ventana de Hann
[S, F, T_spec] = spectrogram(signal, win, win_len - hop, win_len, fs);

%% ── 4. FIGURA 1: Señal temporal + Espectro de magnitud + Fase ───────────────
figure('Name','Transformada de Fourier', 'NumberTitle','off', ...
       'Position',[100 100 1100 700]);

% --- 4a. Señal en el tiempo
subplot(3,2,[1 2]);
plot(t, signal, 'b', 'LineWidth', 1.2);
title('Señal en el dominio temporal');
xlabel('Tiempo (s)');  ylabel('Amplitud');
xlim([0 0.1]);         % primeros 100 ms para mejor visualización
grid on;

% --- 4b. Espectro de magnitud (escala lineal)
subplot(3,2,3);
stem(freqs, magnitude, 'r', 'MarkerSize', 3, 'LineWidth', 0.8);
title('Espectro de magnitud (lineal)');
xlabel('Frecuencia (Hz)');  ylabel('|X(f)|');
xlim([0 500]);
grid on;

% --- 4c. Espectro de magnitud (escala logarítmica dB)
subplot(3,2,4);
plot(freqs, 20*log10(magnitude + 1e-12), 'm', 'LineWidth', 1.2);
title('Espectro de magnitud (dB)');
xlabel('Frecuencia (Hz)');  ylabel('Magnitud (dB)');
xlim([0 500]);  ylim([-80 10]);
grid on;

% --- 4d. Espectro de fase
subplot(3,2,5);
plot(freqs, rad2deg(phase), 'g', 'LineWidth', 1.0);
title('Espectro de fase');
xlabel('Frecuencia (Hz)');  ylabel('Fase (°)');
xlim([0 500]);
grid on;

% --- 4e. Espectrograma 2D (tiempo × frecuencia)
subplot(3,2,6);
imagesc(T_spec, F, 20*log10(abs(S) + 1e-12));
axis xy;
colormap(gca, 'jet');
colorbar;
title('Espectrograma 2D (STFT)');
xlabel('Tiempo (s)');  ylabel('Frecuencia (Hz)');
ylim([0 500]);
clim([-80 0]);

sgtitle('Análisis de Fourier — fs = 1000 Hz', 'FontSize', 14, 'FontWeight', 'bold');

%% ── 5. FIGURA 2: Representación 2D de la energía espectral ──────────────────
figure('Name','Densidad espectral de potencia 2D', 'NumberTitle','off', ...
       'Position',[150 150 900 500]);

% PSD estimada con el método de Welch
[pxx, fpxx] = pwelch(signal, hann(256), 128, 512, fs);

subplot(1,2,1);
area(fpxx, 10*log10(pxx), 'FaceColor',[0.2 0.5 0.8], 'FaceAlpha',0.6, 'EdgeColor','b');
title('Densidad Espectral de Potencia (Welch)');
xlabel('Frecuencia (Hz)');  ylabel('PSD (dB/Hz)');
xlim([0 500]);
grid on;

% Mapa de calor 2D: ventanas × bins de frecuencia (magnitud del espectrograma)
subplot(1,2,2);
pcolor(T_spec, F(F<=500), abs(S(F<=500,:)));
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

threshold = 0.05;
dom_idx   = find(magnitude > threshold);
[~, ord]  = sort(magnitude(dom_idx), 'descend');
dom_idx   = dom_idx(ord);

for k = 1:min(10, length(dom_idx))
    fprintf('%-20.1f %-12.4f\n', freqs(dom_idx(k)), magnitude(dom_idx(k)));
end
