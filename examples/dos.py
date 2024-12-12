#!/usr/bin/env python3

import ebmb
import matplotlib.pyplot as plt

e, dos = ebmb.square_dos('dos.in')

results = ebmb.get(
   tell=False,
   normal=False,
   dos='dos.in',
   mu=0.5 * e[-1],
   lower=1.6 * e[0],
   upper=0.6 * e[-1],
   eta=1e-4,
   resolution=10001,
   )

print('Integral of DOS: %g' % results['states'])

plt.plot(e, dos)
plt.plot(results['omega'] + results['mu'], results['DOS'])

plt.show()
