#!/usr/bin/env python

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

results = ebmb.get(T=0.5 * ebmb.get(program='critical', **parameters),
    **parameters)

fig, ax = plt.subplots(3, sharex='col', sharey='row')

for x, X in zip(ax, ['Z', 'chi', 'Delta']):
    x.plot(results['iomega'], results[X], 'bo')
    x.plot(results['omega'], results['Re[%s]' % X], 'r-')
    x.plot(results['omega'], results['Im[%s]' % X], 'r--')

    x.set_ylabel('$Z$' if X == 'Z' else r'$\%s$ (eV)' % X)
    x.grid()

plt.xlabel(r'$\omega$ (eV)')
plt.xlim(0.0, parameters['upper'])

plt.show()
