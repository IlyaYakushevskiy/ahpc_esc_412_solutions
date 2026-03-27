import numpy as np
import matplotlib.pyplot as plt


files = {
    r"$100^3$ Box": "power_log_PCS_100.txt",
    r"$200^3$ Box": "power_log_PCS_200.txt",
    r"$500^3$ Box": "power_log_PCS_500.txt"
}

plt.figure(figsize=(10, 6))

for label, filename in files.items():
    try:
        k, pk = np.loadtxt(filename, unpack=True)
        plt.loglog(k, pk, marker='o', markersize=4, linestyle='-', alpha=0.8, label=label)
    except FileNotFoundError:
        print(f"Warning: {filename} not found. Skipping.")

plt.title("Logarithmic Power Spectrum $P(k)$")
plt.xlabel("Wavenumber $k$")
plt.ylabel("Power $P(k)$")
plt.grid(True, which="both", ls="--", alpha=0.4)
plt.legend()

plt.tight_layout()
plt.savefig("power_spectrum_log.png", dpi=300)
plt.show()