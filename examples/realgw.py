#!/usr/bin/env python3

import ebmb
import matplotlib.pyplot as plt

dos = 'dos.in'
a2f = 'a2f.in'

ebmb.chain_dos(dos, de=5e-3, t=1.0)
ebmb.chain_a2F(a2f, dw=1e-2, wlog=2.0, l=1.0)

fig, ax = plt.subplots(3, sharex='col')

for (realgw, eta0Im, style, label) in [
    (True, True, 'm', 'real axis (KK)'),
    (True, False, 'c--', 'real axis'),
    (False, False, 'k:', 'imag. axis'),
]:
    results = ebmb.get(
        normal=True,
        Sigma=True,
        chiC=True,
        realgw=realgw,
        eta0Im=eta0Im,
        dos=dos,
        a2F=a2f,
        muC=1.0,
        divdos=True,
        cutoff=10.0,
        lower=-20.0,
        upper=+20.0,
        resolution=501,
        eta=0.1,
        n=1.5,
        T=300.0,
    )

    for x, X in zip(ax, ['Sigma', 'Z', 'chi']):
        x.plot(results['omega'], results['Re[%s]' % X], style, label=label)
        x.plot(results['omega'], results['Im[%s]' % X], style)

ax[0].legend()

ax[0].set_ylabel(r'$\Sigma$ (eV)')
ax[1].set_ylabel(r'$Z$')
ax[2].set_ylabel(r'$\chi$ (eV)')

plt.xlabel(r'$\omega$ (eV)')

plt.savefig('realgw.png')
