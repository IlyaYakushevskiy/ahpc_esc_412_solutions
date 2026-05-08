import matplotlib.pyplot as plt

n = [1000000, 8000000, 125000000]
read = [0.0285946, 0.186947, 3.21436]
assign = [0.00291283, 0.0428932, 1.07015]
proj = [0.000845333, 0.00829096, 0.145411]

fig, axs = plt.subplots(1, 3, figsize=(15, 5))

axs[0].plot(n, read, 'o-')
axs[0].plot(n, assign, 's-')
axs[0].plot(n, proj, '^-')
axs[0].set_title('Linear')

axs[1].plot(n, read, 'o-')
axs[1].plot(n, assign, 's-')
axs[1].plot(n, proj, '^-')
axs[1].set_yscale('log')
axs[1].set_title('Log-Linear')

axs[2].plot(n, read, 'o-', label='Reading')
axs[2].plot(n, assign, 's-', label='Assignment')
axs[2].plot(n, proj, '^-', label='Projection')
axs[2].set_xscale('log')
axs[2].set_yscale('log')
axs[2].set_title('Log-Log')
axs[2].legend()

plt.tight_layout()
plt.savefig("time-plots-task3-2.png")