#!/usr/bin/env python

from numpy import fromfile, int32, float64
from os.path import abspath, dirname, join
from subprocess import call

def load(filename):
    data = {}

    dtype = float64
    shape = [1]

    with open(filename, 'rb') as file:
        while True:
            name = ''.join(iter(lambda: file.read(1) or ':', ':'))

            if name == 'REAL':
                dtype = float64

            elif name == 'INT':
                dtype = int32

            elif name == 'DIM':
                shape = fromfile(file, int32, *fromfile(file, int32, 1))

            elif name:
                data[name] = fromfile(file, dtype,
                    shape.prod()).reshape(shape, order='Fortran')
            else:
                return data

def run(executable=join(dirname(abspath(__file__)), 'eb'),
        filename='~temporary.in', **parameters):

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

            ('epsilon', 1e-15)]: # negligible float difference

            print >> file, parameters.get(parameter, default)

    call([executable, filename])

    return load(filename.rsplit('.', 1)[0] + '.dat')

def squareDOS(name='dos.in', t=0.25, eF=0.5, n=401):
    from scipy import linspace, pi
    from scipy.special import ellipk

    if not n % 2:
        n += 1

    e, de = linspace(-4 * t, 4 * t, n, retstep=True)

    dos = ellipk(1 - (0.25 * e / t) ** 2) / (2 * pi ** 2 * t)

    dos[n // 2] = 0.0
    dos[n // 2] = 1 / de - (dos[0] + 2 * sum(dos[1:-1]) + dos[-1]) / 2

    e += 4 * t - eF

    with open(name, 'w') as file:
        file.write('%d\n\n' % n)

        for i in range(n):
            file.write('% .10f %.10f\n' % (e[i], dos[i]))

if __name__ == '__main__':
    squareDOS()
    print run()
