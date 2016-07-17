#!/usr/bin/env python

from os import path
import numpy as np
import subprocess

try:
    from scipy.special import ellipk
except ImportError:
    print 'Warning: Elliptic integral not available'

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
                    shape.prod()).reshape(shape)
            else:
                return data

def escape(value):
    string = str(value)

    if any(character in string for character in '/ '):
        string = "'%s'" % string.replace("'", "''")

    return string

def new(filename, lamda=1.0, **parameters):
    with open(filename, 'w') as file:
        for parameter, default in [
            ('T', 10.0), # temperature (K)

            ('error', 1e-03), # valid error of critical temperature (K)
            ('bound', 1e+00), # lower bound of critical temperature (K)
            ('small', 1e-10), # maximum gap at critical temperature (eV)

            ('bands', 1), # number of electronic bands

            ('omegaE', 0.02), # Einstein frequency (eV)
            ('lambda', lamda), # electron-phonon coupling
            ('muStar', 0.15), # Coulomb pseudo-potential

            ('DOSfile', 'none'), # file with density of states

            ('upper', 10.0), # overall cutoff (Einstein frequency)
            ('lower', -1.0), # Coulomb cutoff (Einstein frequency)

            ('limit', 100000), # maximum number of fixed-point steps

            ('measurable', True), # find measurable gap?
            ('resolution',  300), # real axis resolution

            ('form', 'data'), # output format
            ('edit', 'ES15.6E3'), # number format

            ('standalone', True), # include parameters in output file?
            ('rescale', True), # rescale Coulomb pseudo-potential?

            ('epsilon', 1e-15)]: # negligible float difference

            value = parameters.get(parameter, default)

            print >> file, '\n'.join(' '.join(map(escape, row))
                for row in value) if np.shape(value) else escape(value)

def run(filename, executable='eb_local'):
    subprocess.call([executable, filename])

def get(infile='~temporary.in', executable='eb_local', replace=True,
    **parameters):

    outfile = infile.rsplit('.', 1)[0] + '.dat'

    if replace or not path.exists(outfile):
        new(infile, **parameters)
        run(infile, executable)

    if path.exists(outfile):
        return load(outfile)

def squareDOS(e, t=0.25):
    return ellipk(1 - (0.25 * e / t) ** 2) / (2 * np.pi ** 2 * t)

def squareDOSGauss(e, t=0.25, n=200, sigma=0.02, cutoff=None):
    k = np.linspace(-np.pi, np.pi, n, endpoint=False)

    e1 = -2 * t * np.cos(k)
    e2 = np.array([I + II for I in e1 for II in e1])

    r = abs(e2 - e)

    if cutoff is not None:
        r = r[np.where(r < cutoff)]

    return sum(np.exp(-(r / sigma) ** 2) / (np.sqrt(np.pi) * sigma)) / n ** 2

def squareDOSfile(name='dos.in', eF=0.5, n=401, t=0.25, bands=1, replace=True):
    if not replace and path.exists(name):
        return

    if not n % 2:
        n += 1

    e, de = np.linspace(-4 * t, 4 * t, n, retstep=True)

    dos = squareDOS(e, t)

    dos[n // 2] = 0.0
    dos[n // 2] = 1 / de - (dos[0] + 2 * sum(dos[1:-1]) + dos[-1]) / 2

    e -= e[0] + eF

    with open(name, 'w') as file:
        file.write('%d\n\n' % n)

        for i in range(n):
            print >> file, '% .10f' % e[i] + ' %.10f' % dos[i] * bands

if __name__ == '__main__':
    np.set_printoptions(threshold=9, edgeitems=1)

    for item in sorted(get().items()):
        print ('%9s = %s' % item).replace('\n', '\n' + ' ' * 12)
