import numpy as np
import matplotlib.pyplot as plt

data = np.fromfile("projection.dat", dtype=np.float32).reshape(100, 100)
plt.imshow(np.log10(data + 1e-10), cmap='magma')
plt.savefig("proj.png")