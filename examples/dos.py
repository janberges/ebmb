#!/usr/bin/env python3

import ebmb
import matplotlib.pyplot as plt

e, dos = ebmb.square_dos('dos.in')

results = ebmb.get(
   tell=False,
   normal=False,
   lamda=2.0,
   dos='dos.in',
   mu=0.5,
   lower=-0.25,
   upper=+0.25,
   eta=1e-4,
   resolution=10001,
   stable=True,
   )

print('Integral of DOS: %g' % results['states'])

plt.fill_between(e, dos, color='lightgray')
plt.plot(results['omega'] + results['mu'], results['DOS'])

plt.show()
