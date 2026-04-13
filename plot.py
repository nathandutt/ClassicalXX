import csv
import numpy as np
import matplotlib.pyplot as plt
import sys

def read_config(filename):
    config = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) == 2:
                key, value = parts
                try:
                    config[key] = float(value)
                except ValueError:
                    print(f"Warning: Skipping invalid line (non-numeric value) in config: {line}")
            else:
                print(f"Warning: Skipping invalid line (expected key-value pair) in config: {line}")
    return config
config = read_config('config.txt')
x_range = config.get('x_range', 20)
theta_l = config.get('theta_l', 1.)
theta_r = config.get('theta_r', 1.)
S_z_r = np.cos(np.pi*theta_r)
S_z_l = np.cos(np.pi*theta_l)
dx = config.get('dx', 0.1)
x = np.arange(-x_range, x_range, dx)
x = x[:int(2 * x_range / dx)]

fig, ax = plt.subplots(figsize=(10,10))
line1, = ax.plot(x, np.zeros_like(x), label=r'$\rho$', color = "black")
line2, = ax.plot(x, np.zeros_like(x), label=r'$\phi_x$', color = "green")
time_text = ax.text(0.05, 0.05, '', transform=ax.transAxes, fontsize=16)

ax.set_xlabel('x', fontsize=20)
ax.set_ylabel(r'$S_z$', fontsize = 20)
plt.xlim(-x_range, x_range)
ax.set_xticks([-x_range*0.9, 0, x_range*0.9])
ax.set_yticks([-1, -0.5, 0., 0.5, 1.])
plt.ylim(-1.1, 1.1)
ax.grid(False)
plt.title(rf'Plot of $S_z$ (Initially $S_{{\mathrm{{z, right}}}} = {S_z_r:.2f}$, $S_{{\mathrm{{z, left}}}} = {S_z_l:.2f}$)')
plt.ion()
with open('evolution.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    row_iter = iter(reader)
    for row_1 in row_iter:
        row_2 = next(row_iter, None)
        row_3 = next(row_iter, None)

        S_x = [float(x) for x in row_1]
        S_y = [float(x) for x in row_2]
        S_z = [float(x) for x in row_3]

        t = S_z[0]
        S_x = np.array(S_x[1:])
        S_y=np.array(S_y[1:])
        S_z = np.array(S_z[1:])

        compl = (S_x + 1j * S_y)/(np.sqrt((1-(S_z**2)))+1e-8)
        compl[S_z > 0.999] = 1.
        compl[S_z < -0.999] = 1.
        
        v = np.imag(np.gradient(compl, dx)/compl)
        theta = np.arccos(S_z)/np.pi
        sigma = np.arccos(v)/np.pi
        rplus = theta+sigma
        rminus = theta - sigma
        line1.set_ydata(S_z)
        line2.set_ydata(v)
        time_text.set_text(rf'$\tau$={t:.2f}')
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.waitforbuttonpress()


plt.ioff()
plt.show()

