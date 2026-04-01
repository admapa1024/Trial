import numpy as np
import matplotlib.pyplot as plt


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


def plot_results(t, signal, freqs, magnitude):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

    ax1.plot(t[:200], signal[:200])
    ax1.set_title("Señal en el tiempo")
    ax1.set_xlabel("Tiempo (s)")
    ax1.set_ylabel("Amplitud")
    ax1.grid(True)

    ax2.plot(freqs, magnitude)
    ax2.set_title("Espectro de frecuencias (FFT)")
    ax2.set_xlabel("Frecuencia (Hz)")
    ax2.set_ylabel("Magnitud")
    ax2.set_xlim(0, 500)
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig("fourier_result.png", dpi=150)
    plt.show()
    print("Gráfica guardada en fourier_result.png")


def print_peaks(freqs, magnitude, threshold=0.1):
    """Imprime las frecuencias dominantes."""
    peaks = [(freqs[i], magnitude[i]) for i in range(len(magnitude)) if magnitude[i] > threshold]
    peaks.sort(key=lambda x: -x[1])
    print("\nFrecuencias dominantes:")
    print(f"{'Frecuencia (Hz)':>18}  {'Magnitud':>10}")
    print("-" * 32)
    for freq, mag in peaks[:10]:
        print(f"{freq:>18.1f}  {mag:>10.4f}")


if __name__ == "__main__":
    FS = 1000       # frecuencia de muestreo (Hz)
    DURATION = 1.0  # duración de la señal (segundos)

    t, signal = generate_signal(fs=FS, duration=DURATION)
    freqs, magnitude = compute_fft(signal, fs=FS)

    print_peaks(freqs, magnitude)
    plot_results(t, signal, freqs, magnitude)
