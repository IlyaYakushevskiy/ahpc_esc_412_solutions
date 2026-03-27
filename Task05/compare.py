import numpy as np
import matplotlib.pyplot as plt
import pathlib

BASE = pathlib.Path(__file__).resolve().parent

left_compare = BASE / "../Task03/ex03(ex02cp)/projection.dat"
right_compare = BASE / "density_NGP_100.bin"

projected_left = np.fromfile(left_compare, dtype=np.float32)
nGrid_left = int(np.sqrt(projected_left.size))
projected_left = projected_left.reshape((nGrid_left, nGrid_left))

with open(right_compare, "rb") as f:
    nGrid_right = np.fromfile(f, dtype=np.int32, count=1)[0]
    projected_right = np.fromfile(
        f, dtype=np.float32, count=nGrid_right * nGrid_right
    ).reshape((nGrid_right, nGrid_right))

factor = nGrid_left // nGrid_right
projected_left = projected_left.reshape(
    nGrid_right, factor, nGrid_right, factor
).mean(axis=(1, 3))

errors = np.abs(projected_left - projected_right)
print("Max error:", errors.max())

plt.imshow(np.log10(projected_left + 1e-10), cmap="magma")
plt.savefig("proj_left.png")
plt.close()

plt.imshow(np.log10(projected_right + 1e-10), cmap="magma")
plt.savefig("proj_right.png")
plt.close()