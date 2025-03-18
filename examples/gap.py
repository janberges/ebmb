#!/usr/bin/env python3

import ebmb
import matplotlib.pyplot as plt
import numpy as np

para = dict(tell=False, lamda=[[1, 1e-3], [1e-3, 2]])

Tc = ebmb.get('critical', **para)
T = np.linspace(0.1, 1.1, 100) * Tc

Delta = np.empty((T.size, 2))

for i, para['T'] in enumerate(T):
    print('T = %(T)g K' % para)

    Delta[i] = ebmb.get(**para)['Delta'][:, 0]

plt.plot(T, 1e3 * Delta)

plt.xlabel(r'$T$ (K)')
plt.ylabel(r'$\Delta(\mathrm{i} \omega_0)$ (meV)')
plt.grid()

plt.savefig('gap.png')
