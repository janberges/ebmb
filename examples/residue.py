#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-10.0, 10.0, 101)
kT = 1.0

fig, sub = plt.subplots(3, 3)

for N, ax in enumerate(sub.flat):
    summation = -kT * sum(1 / (1j * (2 * n + 1) * np.pi * kT - x)
        for n in range(-N, N)).real

    residue = np.arctan2(x, (2 * N + 1) * np.pi * kT) / np.pi

    ax.fill_between(x, summation, summation + residue)
    ax.plot(x, np.tanh(x / (2 * kT)) / 2, 'k')

    ax.label_outer()

plt.show()
