#!/usr/bin/env python3

import ebmb
import matplotlib.pyplot as plt

dos = 'dos.in'
a2f = 'a2f.in'

ebmb.chain_dos(dos, de=5e-3, t=1.0)
ebmb.chain_a2F(a2f, dw=1e-2, l=1.0, wlog=2.0)

fig, ax = plt.subplots(3, 2, sharex='col', sharey='row')

for (realgw, krakro, style, label) in [
    (True, True, 'm', 'Kramers-Kronig'),
    (True, False, 'c--', 'direct'),
    (False, False, 'k:', 'Pad\xe9'),
]:
    results = ebmb.get(
        normal=True,
        chiC=True,
        realgw=realgw,
        krakro=krakro,
        dos=dos,
        a2F=a2f,
        muC=1.0,
        cutoff=10.0,
        lower=-20.0,
        upper=+20.0,
        resolution=501,
        eta=0.1,
        n=1.5,
        T=300.0,
    )

    for x, X in zip(ax, ['Sigma', 'Z', 'chi']):
        x[0].plot(results['omega'], results['Re[%s]' % X], style, label=label)
        x[0].plot(results['omega'], results['Im[%s]' % X], style)

        if X == 'Sigma':
            x[1].plot(results['iomega'], results['chi'], style, label=label)
            x[1].plot(results['iomega'], results['domega'], style)
        else:
            x[1].plot(results['iomega'], results[X], style, label=label)

        x[0].set_ylabel('$Z$' if X == 'Z' else r'$\%s$ (eV)' % X)

        x[0].grid()
        x[1].grid()

ax[-1, 0].set_xlabel(r'$\omega$ (eV)')
ax[-1, 1].set_xlabel(r'$\omega_n$ (eV)')

ax[1, 1].legend()

plt.savefig('realgw.png')
