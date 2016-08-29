#!/usr/bin/env python

"""Wrapper and auxiliary functions for Eliashberg solver ebmb"""

import itertools
import numpy as np
from os import path
import subprocess

try:
    from scipy.special import ellipk
except ImportError:
    print 'square_dos not available'

def get(program='ebmb', file='~temporary.dat', replace=True, **parameters):
    """Run 'ebmb', 'tc' or 'critical' and load results.

    Parameters
    ----------
    program : str
        Name of or path to executable.
    file : str
        Path to output file.
    replace : bool
        Overwrite existing output file?
    **parameters
        Program parameters.

    Returns
    -------
    dict
        Returned if `program` corresponds to 'ebmb'.
        Self-energy components etc.
    ndarray
        Returned otherwise.
        Critical parameter(s).
    """
    if replace or not path.exists(file):
        run(program, file=file, **parameters)

    if program.endswith('ebmb'):
        return load(file)
    else:
        return load_floats(file)

def run(program='ebmb', **parameters):
    """Run 'ebmb', 'tc' or 'critical'.

    Parameters
    ----------
    program : str
        Name of or path to executable.
    **parameters
        Program parameters.
    """
    command = [program]

    for key, value in parameters.items():
        command.append('='.join([key, ','.join(map(str, np.ravel(value)))]))

    subprocess.call(command)

def load(filename):
    """Load output file of 'ebmb'.

    Parameters
    ----------
    filename : str
        Path to output file.

    Returns
    -------
    dict
        Self-energy components etc.
    """
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

def load_floats(filename):
    """Load output file of 'tc' or 'critical'.

    Parameters
    ----------
    filename : str
        Path to output file.

    Returns
    -------
    ndarray
        Critical parameter(s).
    """
    with open(filename, 'rb') as file:
        data = np.fromfile(file, np.float64)

    return data if data.size > 1 else data[0]

def dos(file, epsilon, domain, filters=[], n=101, replace=True):
    """Calculate subdomain-resolved density of states and save it to file.

    Parameters
    ----------
    file : str
        Path to output file.
    epsilon : function
        Band structure.
    domain : list of ndarray
        Discretized domains of arguments of `epsilon`.
    filters : list of function
        N filters defining N + 1 subdomains.
    n : int
        Resolution of density of states.
    replace : bool
        Overwrite existing output file?

    Returns
    -------
    ndarray
        Energy.
    ndarray
        Subdomain-resolved density of states.
    """
    if not replace and path.exists(file):
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

    return e, dos if p > 1 else dos[:, 0]

def square_dos(file='dos.in', n=401, t=0.25, replace=True):
    """Calculate density of states of square lattice and save it to file.

    Parameters
    ----------
    file : str
        Path to output file.
    n : int
        Resolution of density of states.
    t : float
        Hopping parameter.
    replace : bool
        Overwrite existing output file?

    Returns
    -------
    ndarray
        Energy.
    ndarray
        Density of states.
    """
    if not replace and path.exists(file):
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

    return e, dos

if __name__ == '__main__':
    np.set_printoptions(threshold=9, edgeitems=1)

    square_dos('dos.in')

    for item in sorted(get(dos='dos.in', n=0.5, tell=False).items()):
        print ('%9s = %s' % item).replace('\n', '\n' + ' ' * 12)
