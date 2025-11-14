#!/usr/bin/env python3

import ebmb
import matplotlib.pyplot as plt

parameters = dict(
    tell=False,
    lamda=1.0,
    omegaE=0.025,
    dos='dos.in',
    n=0.5,
    cutoff=100.0,
    upper=0.2,
    points=1000,
)

ebmb.square_dos(parameters['dos'])

Tc = ebmb.get('critical', **parameters)

results = ebmb.get(T=0.5 * Tc, **parameters)

fig, ax = plt.subplots(3, 2, sharex='col', sharey='row')

for x, X in zip(ax, ['Z', 'chi', 'Delta']):
    x[0].plot(results['omega'], results['Re[%s]' % X], 'r-')
    x[0].plot(results['omega'], results['Im[%s]' % X], 'r--')
    x[1].plot(results['iomega'], results[X], 'bo')

    x[0].set_ylabel('$Z$' if X == 'Z' else r'$\%s$ (eV)' % X)
    x[1].set_xlim(0.0, parameters['upper'])

    x[0].grid()
    x[1].grid()

ax[-1, 0].set_xlabel(r'$\omega$ (eV)')
ax[-1, 1].set_xlabel(r'$\omega_n$ (eV)')

fig.savefig('selfen.png')
