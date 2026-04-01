import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys


def load_audio_file(path):
    """Carga un archivo de audio y devuelve t, señal y fs.

    Formatos soportados:
      - soundfile (pip install soundfile): WAV, FLAC, OGG, AIFF, AU, etc.
      - scipy (pip install scipy):         WAV únicamente (fallback).
      - MP3: requiere instalar pydub y ffmpeg (pip install pydub).
    """
    ext = path.rsplit(".", 1)[-1].lower()

    # ── MP3: requiere pydub + ffmpeg ─────────────────────────────────────────
    if ext == "mp3":
        try:
            from pydub import AudioSegment
        except ImportError:
            print("Error: para leer MP3 instala pydub:  pip install pydub")
            print("       y asegúrate de tener ffmpeg en el PATH.")
            sys.exit(1)
        audio = AudioSegment.from_mp3(path)
        fs = audio.frame_rate
        samples = np.array(audio.get_array_of_samples(), dtype=np.float64)
        if audio.channels > 1:
            samples = samples.reshape(-1, audio.channels).mean(axis=1)
        samples /= 2 ** (audio.sample_width * 8 - 1)
        t = np.linspace(0, len(samples) / fs, len(samples), endpoint=False)
        return t, samples, fs

    # ── Resto de formatos: soundfile (WAV, FLAC, OGG, AIFF, AU…) ────────────
    try:
        import soundfile as sf
        data, fs = sf.read(path, dtype="float64", always_2d=False)
        if data.ndim > 1:
            data = data.mean(axis=1)
        t = np.linspace(0, len(data) / fs, len(data), endpoint=False)
        return t, data, fs
    except ImportError:
        pass  # soundfile no disponible, intentar con scipy

    # ── Fallback: scipy (solo WAV) ───────────────────────────────────────────
    try:
        from scipy.io import wavfile
    except ImportError:
        print("Error: instala soundfile (recomendado) o scipy:")
        print("       pip install soundfile")
        sys.exit(1)

    fs, data = wavfile.read(path)
    if data.dtype == np.int16:
        data = data.astype(np.float64) / 32768.0
    elif data.dtype == np.int32:
        data = data.astype(np.float64) / 2147483648.0
    elif data.dtype == np.uint8:
        data = (data.astype(np.float64) - 128.0) / 128.0
    else:
        data = data.astype(np.float64)
    if data.ndim > 1:
        data = data.mean(axis=1)
    t = np.linspace(0, len(data) / fs, len(data), endpoint=False)
    return t, data, fs


def generate_signal(fs=1000, duration=1.0):
    """Genera una señal de ejemplo compuesta por varias frecuencias."""
    t = np.linspace(0, duration, int(fs * duration), endpoint=False)
    # Señal = 3 componentes sinusoidales: 50 Hz, 120 Hz y 300 Hz
    signal = (
        1.0 * np.sin(2 * np.pi * 50 * t) +
        0.5 * np.sin(2 * np.pi * 120 * t) +
        0.25 * np.sin(2 * np.pi * 300 * t)
    )
    return t, signal


def compute_fft(signal, fs):
    """Calcula la FFT y devuelve frecuencias y magnitudes."""
    n = len(signal)
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    fft_vals = np.fft.rfft(signal)
    magnitude = np.abs(fft_vals) / n * 2  # normalización (x2 por espectro unilateral)
    return freqs, magnitude


def plot_results(t, signal, freqs, magnitude, fs):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

    # Mostrar los primeros 100 ms de la señal
    samples_100ms = min(len(t), int(fs * 0.1))
    ax1.plot(t[:samples_100ms], signal[:samples_100ms])
    ax1.set_title("Señal en el tiempo (primeros 100 ms)")
    ax1.set_xlabel("Tiempo (s)")
    ax1.set_ylabel("Amplitud")
    ax1.grid(True)

    ax2.plot(freqs, magnitude)
    ax2.set_title("Espectro de frecuencias (FFT)")
    ax2.set_xlabel("Frecuencia (Hz)")
    ax2.set_ylabel("Magnitud")
    ax2.set_xlim(0, fs / 2)
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig("fourier_result.png", dpi=150)
    plt.show()
    print("Gráfica guardada en fourier_result.png")


def print_peaks(freqs, magnitude, threshold=None):
    """Imprime las frecuencias dominantes."""
    if threshold is None:
        threshold = max(magnitude) * 0.05  # 5% del pico máximo
    peaks = [(freqs[i], magnitude[i]) for i in range(len(magnitude)) if magnitude[i] > threshold]
    peaks.sort(key=lambda x: -x[1])
    print("\nFrecuencias dominantes:")
    print(f"{'Frecuencia (Hz)':>18}  {'Magnitud':>10}")
    print("-" * 32)
    for freq, mag in peaks[:10]:
        print(f"{freq:>18.1f}  {mag:>10.4f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Análisis de Fourier de una señal.",
        epilog="Si no se proporciona archivo, se usa una señal sintética de ejemplo."
    )
    parser.add_argument(
        "audio_file",
        nargs="?",
        default=None,
        help="Ruta al archivo de audio (WAV, FLAC, OGG, AIFF, MP3…)"
    )
    args = parser.parse_args()

    if args.audio_file:
        print(f"Cargando archivo de audio: {args.audio_file}")
        t, signal, FS = load_audio_file(args.audio_file)
        print(f"  Frecuencia de muestreo : {FS} Hz")
        print(f"  Duración               : {len(signal)/FS:.2f} s")
        print(f"  Muestras               : {len(signal)}")
    else:
        print("No se especificó archivo de audio. Usando señal sintética de ejemplo.")
        FS = 1000
        DURATION = 1.0
        t, signal = generate_signal(fs=FS, duration=DURATION)

    freqs, magnitude = compute_fft(signal, fs=FS)
    print_peaks(freqs, magnitude)
    plot_results(t, signal, freqs, magnitude, FS)
