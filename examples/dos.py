#!/usr/bin/env python3

import ebmb
import matplotlib.pyplot as plt

e, dos = ebmb.square_dos('dos.in')

mu = 0.5

results = ebmb.get(
   tell=False,
   normal=False,
   lamda=2.0,
   dos='dos.in',
   mu=mu,
   lower=-0.25,
   upper=+0.25,
   eta=1e-4,
   points=10001,
   stable=True,
)

print('Integral of noninteracting DOS: %g' % results['states'])
print('Integral of quasi-particle DOS: %g' % results['inspect'])

plt.fill_between(e - mu, dos, color='lightgray')
plt.plot(results['omega'], results['DOS'])

plt.savefig('dos.png')
