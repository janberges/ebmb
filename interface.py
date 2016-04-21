#!/usr/bin/env python

import numpy as np
import subprocess

try:
    from scipy.special import ellipk
except ImportError:
    print 'Warning: squareDOS not available'

def load(filename):
    data = {}

    with open(filename, 'rb') as file:
        while True:
            name = ''.join(iter(lambda: file.read(1) or ':', ':'))

            if name == 'REAL':
                dtype = np.float64

            elif name == 'INT':
                dtype = np.int32

            elif name == 'DIM':
                shape = np.fromfile(file, np.int32,
                    *np.fromfile(file, np.int32, 1))

            elif name:
                data[name] = np.fromfile(file, dtype,
                    shape.prod()).reshape(shape, order='Fortran')
            else:
                return data

def run(executable='eb_local', filename='~temporary.in', **parameters):

    if 'DOSfile' in parameters:
        if any(c in parameters['DOSfile'] for c in '/ '):
            parameters['DOSfile'] = "'%s'" % parameters['DOSfile'].replace(
                "'", "''")

    with open(filename, 'w') as file:

        for parameter, default in [
            ('T', 10.0), # temperature (K)

            ('small', 1e-10), # negligible gap (eV)
            ('error', 1e-10), # error of critical temperature (K)

            ('omegaE', 0.020), # Einstein frequency (eV)
            ('lambda', 1.748), # electron-phonon coupling
            ('muStar', 0.100), # Coulomb pseudo-potential

            ('DOSfile', 'none'), # file with density of states

            ('upper', 10.0), # overall cutoff frequency (eV)
            ('lower',  5.0), # Coulomb cutoff frequency (eV)

            ('limit', 100000), # maximum number of fixed-point steps

            ('measurable', True), # find measurable gap?
            ('resolution',  300), # real axis resolution

            ('form', 'data'), # output format
            ('standalone', True), # include parameters in output file?

            ('epsilon', 1e-15)]: # negligible float difference

            print >> file, parameters.get(parameter, default)

    subprocess.call([executable, filename])

    return load(filename.rsplit('.', 1)[0] + '.dat')

def squareDOS(name='dos.in', t=0.25, eF=0.5, n=401):
    if not n % 2:
        n += 1

    e, de = np.linspace(-4 * t, 4 * t, n, retstep=True)

    dos = ellipk(1 - (0.25 * e / t) ** 2) / (2 * np.pi ** 2 * t)

    dos[n // 2] = 0.0
    dos[n // 2] = 1 / de - (dos[0] + 2 * sum(dos[1:-1]) + dos[-1]) / 2

    e += 4 * t - eF

    with open(name, 'w') as file:
        file.write('%d\n\n' % n)

        for i in range(n):
            file.write('% .10f %.10f\n' % (e[i], dos[i]))

if __name__ == '__main__':
    np.set_printoptions(threshold=3, edgeitems=1)

    for item in sorted(run().items()):
        print '%14s = %s' % item
