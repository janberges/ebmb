#!/usr/bin/env python

import itertools
import numpy as np
from os import path
import subprocess

try:
    from scipy.special import ellipk
except ImportError:
    print "'squareDOSfile' not available"

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

def run(program='ebmb', **parameters):
    command = [program]

    for key, value in parameters.items():
        command.append('='.join([key, ','.join(map(escape, np.ravel(value)))]))

    subprocess.call(command)

def get(file='~temporary.dat', replace=True, program='ebmb', **parameters):
    if replace or not path.exists(file):
        run(program, file=file, **parameters)

    if program.endswith('ebmb'):
        return load(file)
    else:
        with open(file, 'rb') as file:
            data = np.fromfile(file, np.float64)

        return data if data.size > 1 else data[0]

def ebmb(**parameters):
    return get(program='ebmb', **parameters)

def tc(**parameters):
    return get(program='tc', **parameters)

def critical(**parameters):
    return get(program='critical', **parameters)

def squareDOSfile(file='dos.in', n=401, t=0.25, replace=True):
    if not replace and path.exists(name):
        return

    if not n % 2:
        n += 1

    e, de = np.linspace(-4 * t, 4 * t, n, retstep=True)

    dos = ellipk(1 - (e / (4 * t)) ** 2) / (2 * np.pi ** 2 * t)

    dos[n // 2] = 0.0
    dos[n // 2] = 1 / de - (dos[0] + 2 * sum(dos[1:-1]) + dos[-1]) / 2

    with open(file, 'w') as out:
        for i in range(n):
            out.write('% .10f %.10f\n' % (e[i], dos[i]))

def DOSfile(file, epsilon, domain, filters=[], n=101, replace=True):
    if not replace and path.exists(name):
        return

    points = np.prod(map(len, domain))

    pocket = np.empty(points, dtype=int)
    energy = np.empty(points)

    for i, x in enumerate(itertools.product(*domain)):
        energy[i] = epsilon(*x)
        pocket[i] = 0

        for element in filters:
            if element(*x): break
            pocket[i] += 1

    emin = energy.min()
    emax = energy.max()

    binned = ((n - 1) * (energy - emin) / (emax - emin)).round().astype(int)

    p = len(filters) + 1

    count = np.zeros((n, p), dtype=int)

    for i in range(points):
        count[binned[i], pocket[i]] += 1

    e, de = np.linspace(emin, emax, n, retstep=True)

    dos = count / (de * count.sum())
    dos[(0, -1), :] *= 2

    with open(file, 'w') as out:
        for i in range(n):
            out.write('% .10f' % e[i])

            for j in range(p):
                out.write(' %.10f' % dos[i, j])

            out.write('\n')

if __name__ == '__main__':
    np.set_printoptions(threshold=9, edgeitems=1)

    for item in sorted(ebmb(tell=False).items()):
        print ('%9s = %s' % item).replace('\n', '\n' + ' ' * 12)
