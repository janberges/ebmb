#!/usr/bin/env python3

import ebmb

lamda = 1.0
omegaE = 0.02
sigma = 0.002

w, a2F = ebmb.gaussian_a2F('a2f.in', l=lamda, w0=omegaE, s=sigma)

Tc_a2F = ebmb.get(
   program='critical',
   tell=False,
   a2F='a2f.in',
)

print('Tc from Eliashberg function: %g K' % Tc_a2F)

Tc_E = ebmb.get(
   program='critical',
   tell=False,
   lamda=lamda,
   omegaE=omegaE,
)

print('Tc from Einstein frequency: %g K' % Tc_E)
