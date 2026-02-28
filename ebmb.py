#!/usr/bin/env python3

# Copyright (C) 2016-2026 Jan Berges
# This program is free software under the terms of the GNU GPLv3 or later.

"""Wrapper and auxiliary functions for Eliashberg solver ebmb"""

__version__ = '2.0.0'

import itertools
import numpy as np
from os import path
import subprocess

try:
    from scipy.special import ellipk
except ImportError:
    print('square_dos not available')

def get(program='ebmb', file=None, replace=True, **parameters):
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
    if file is None:
        file = '~ebmb.tmp.dat'
        shm = '/dev/shm/'

        if path.exists(shm):
            file = shm + file

    if replace or not path.exists(file):
        run(program, file=file, **parameters)

    if program.endswith('ebmb'):
        return load(file)
    else:
        return load_floats(file)

def run(program='ebmb', redirect=False, **parameters):
    """Run 'ebmb', 'tc' or 'critical'.

    Parameters
    ----------
    program : str
        Name of or path to executable. If not found, the path to a directory
        'bin' located next to this module is prepended before trying again.
    redirect : bool
        Do not print but return standard output of program.
    **parameters
        Program parameters.
    """
    command = [program]

    for key, value in parameters.items():
        if value is not None:
            command.append('='.join([key, ','.join(map(str, np.ravel(value)))]))

    for attempt in range(2):
        try:
            if redirect:
                return subprocess.check_output(command)
            else:
                return subprocess.call(command)
        except FileNotFoundError:
            command[0] = path.join(path.dirname(path.abspath(__file__)), 'bin',
                command[0])

def load(file):
    """Load output file of 'ebmb'.

    Parameters
    ----------
    file : str
        Path to output file.

    Returns
    -------
    dict
        Self-energy components etc.
    """
    data = {}

    with open(file, 'rb') as file:
        while True:
            name = ''.join(iter(lambda: chr(ord(file.read(1) or ':')), ':'))

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

def load_floats(file):
    """Load output file of 'tc' or 'critical'.

    Parameters
    ----------
    file : str
        Path to output file.

    Returns
    -------
    ndarray
        Critical parameter(s).
    """
    with open(file, 'rb') as file:
        data = np.fromfile(file, np.float64)

    return data if data.size > 1 else data[0]

def dos(file, epsilon, domain, filters=[], points=101, replace=True):
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
    points : int
        Number of energy points.
    replace : bool
        Overwrite existing output file?

    Returns
    -------
    ndarray
        Energy.
    ndarray
        Subdomain-resolved density of states.

    Examples
    --------
    .. code-block:: python

        e, dos = ebmb.dos('dos.in', epsilon=lambda *k: -np.cos(k).sum() / 2,
            domain=[np.linspace(-np.pi, np.pi, 1000, endpoint=False)] * 2,
            filters=[lambda *k: np.pi ** 2 / 2 <= np.dot(k, k) <= np.pi ** 2])
    """
    if not replace and path.exists(file):
        return

    states = np.prod(list(map(len, domain)))

    energy = np.empty(states)
    pocket = np.empty(states, dtype=int)

    for i, x in enumerate(itertools.product(*domain)):
        energy[i] = epsilon(*x)
        pocket[i] = 0

        for element in filters:
            if element(*x): break
            pocket[i] += 1

    emin = energy.min()
    emax = energy.max()

    binned = ((points - 1)
        * (energy - emin) / (emax - emin)).round().astype(int)

    pockets = len(filters) + 1

    count = np.zeros((points, pockets), dtype=int)

    for i in range(states):
        count[binned[i], pocket[i]] += 1

    e, de = np.linspace(emin, emax, points, retstep=True)

    dos = count / (de * count.sum())
    dos[(0, -1), :] *= 2

    with open(file, 'w') as out:
        for i in range(points):
            out.write('% .10f' % e[i])

            for j in range(pockets):
                out.write(' %.10f' % dos[i, j])

            out.write('\n')

    return e, dos if pockets > 1 else dos[:, 0]

def chain_dos(file='dos.in', de=1e-3, t=0.25, bandwidth=None, replace=True):
    """Calculate density of states of 1D lattice and save it to file.

    Parameters
    ----------
    file : str
        Path to output file.
    de : float
        Energy resolution.
    t : float
        Hopping parameter.
    bandwidth : float
        Alternatively, bandwith.
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

    if bandwidth is not None:
        t = bandwidth / 4

    points = int(round(4 * t / de)) + 1
    points += 1 - points % 2

    e, de = np.linspace(-2 * t, 2 * t, points, retstep=True)

    dos = np.empty(points)

    dos[1:-1] = 1 / (np.pi * np.sqrt((2 * t) ** 2 - e[1:-1] ** 2))

    dos[0] = dos[-1] = 1 / de - sum(dos[1:-1])

    with open(file, 'w') as out:
        for i in range(points):
            out.write('% .10f %.10f\n' % (e[i], dos[i]))

    return e, dos

def chain_a2F(file='a2F.in', dw=1e-4, l=1.0, w0=0.02, wlog=None, replace=True):
    """Calculate and save Eliashberg spectral function of 1D lattice.

    See Eqs. (63) and (64) of Phys. Rev. X 13, 041009 (2023). Averaging the
    electron-phonon coupling over k yields this Eliashberg spectral function.
    The lambda parameter `l` is g0 / w0 squared times the density of states.

    Parameters
    ----------
    file : str
        Path to output file.
    dw : float
        Energy resolution.
    l : float
        Effective electron-phonon coupling.
    w0 : float
        Second-moment average phonon frequency.
    wlog : float
        Logarithmic average phonon frequency.
    replace : bool
        Overwrite existing output file?

    Returns
    -------
    ndarray
        Energy.
    ndarray
        Eliashberg spectral function.
    """
    if not replace and path.exists(file):
        return

    if wlog is None:
        wlog = w0 / np.sqrt(2)

    points = int(round(2 * wlog / dw)) + 1

    w, dw = np.linspace(2 * wlog, 0, points, endpoint=False, retstep=True)
    w = w[::-1]
    dw *= -1

    a2F = np.empty(points)

    a2F[:-1] = l / np.pi * w[:-1] / np.sqrt((2 * wlog) ** 2 - w[:-1] ** 2)

    a2F[-1] = (l / dw - a2F[0] / w[0] - 2 * sum(a2F[1:-1] / w[1:-1])) * w[-1]

    with open(file, 'w') as out:
        for i in range(points):
            out.write('% .10f %.10f\n' % (w[i], a2F[i]))

    return w, a2F

def square_dos(file='dos.in', de=1e-3, t=0.25, bandwidth=None, replace=True):
    """Calculate density of states of square lattice and save it to file.

    Parameters
    ----------
    file : str
        Path to output file.
    de : float
        Energy resolution.
    t : float
        Hopping parameter.
    bandwidth : float
        Alternatively, bandwith.
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

    if bandwidth is not None:
        t = bandwidth / 8

    points = int(round(8 * t / de)) + 1
    points += 1 - points % 2

    e, de = np.linspace(-4 * t, 4 * t, points, retstep=True)

    mid = points // 2

    dos = np.empty(points)

    dos[:mid] = ellipk(1 - (e[:mid] / (4 * t)) ** 2) / (2 * np.pi ** 2 * t)
    dos[mid + 1:] = dos[mid - 1::-1]

    dos[mid] = 0.0
    dos[mid] = 1 / de - dos[0] / 2 - sum(dos[1:-1]) - dos[-1] / 2

    with open(file, 'w') as out:
        for i in range(points):
            out.write('% .10f %.10f\n' % (e[i], dos[i]))

    return e, dos

def box_dos(file='dos.in', de=1e-3, t=0.25, bandwidth=None,
        replace=True):
    """Calculate rectangular density of states and save it to file.

    Parameters
    ----------
    file : str
        Path to output file.
    de : float
        Energy resolution.
    t : float
        One eighth of the bandwidth.
    bandwidth : float
        Alternatively, bandwith.
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

    if bandwidth is not None:
        t = bandwidth / 8

    points = int(round(8 * t / de)) + 1

    e = np.linspace(-4 * t, 4 * t, points)

    dos = np.empty(points)
    dos[:] = 0.125 / t

    with open(file, 'w') as out:
        for i in range(points):
            out.write('% .10f %.10f\n' % (e[i], dos[i]))

    return e, dos

def steplike_dos(file='dos.in', de=1e-3, t=0.25, bandwidth=None, ratio=6.0,
        d=0.02, replace=True):
    """Calculate and save steplike DOS [Akashi, Arita, PRB 88, 014514 (2013)]

    Parameters
    ----------
    file : str
        Path to output file.
    de : float
        Energy resolution.
    t : float
        One eighth of the bandwidth.
    bandwidth : float
        Alternatively, bandwith.
    ratio : float
        Quotient of densities of states after and before the step (N+/N-)
    d : float
        Width of the step.
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

    if bandwidth is not None:
        t = bandwidth / 8

    inner = max(2, int(round(d / de)) + 1)
    outer = max(1, int(round((4 * t - 0.5 * d) / de)))

    points = inner + 2 * outer

    e = np.empty(points)
    dos = np.empty(points)

    e[:+outer] = np.linspace(-4 * t, -0.5 * d, outer, endpoint=False)
    e[-outer:] = np.linspace(4 * t, 0.5 * d, outer, endpoint=False)[::-1]

    e[outer:-outer] = np.linspace(-0.5 * d, 0.5 * d, inner)

    N0 = 0.125 / t
    delta = (ratio - 1.0) / (ratio + 1.0) * N0

    dos[:+outer] = N0 - delta
    dos[-outer:] = N0 + delta

    dos[outer:-outer] = np.linspace(N0 - delta, N0 + delta, inner)

    with open(file, 'w') as out:
        for i in range(points):
            out.write('% .10f %.10f\n' % (e[i], dos[i]))

    return e, dos

def gaussian_a2F(file='a2F.in', dw=1e-4, l=1.0, w0=0.02, s=0.002, wlog=None,
        replace=True):
    """Calculate and save Gaussian example Eliashberg spectral function.

    Parameters
    ----------
    file : str
        Path to output file.
    dw : float
        Energy resolution.
    l : float
        Prefactor (expected value of electron-phonon coupling).
    w0 : float
        Center (expected value of Einstein frequency).
    s : float
        Broadening.
    wlog : float
        Logarithmic average phonon frequency. This is an alias for `w0`.
    replace : bool
        Overwrite existing output file?

    Returns
    -------
    ndarray
        Energy.
    ndarray
        Eliashberg spectral function.
    """
    if not replace and path.exists(file):
        return

    if wlog is not None:
        w0 = wlog

    w = np.arange(w0 - 5 * s, w0 + 5 * s, dw)
    w = w[w >= 0]

    a2F = l * w0 / (2 * np.sqrt(np.pi) * s) * np.exp(-((w - w0) / s) ** 2)

    with open(file, 'w') as out:
        for i in range(len(w)):
            out.write('% .10f %.10f\n' % (w[i], a2F[i]))

    return w, a2F

if __name__ == '__main__':
    np.set_printoptions(threshold=9, edgeitems=1)

    square_dos('dos.in')

    for item in sorted(get(dos='dos.in', n=0.5, tell=False).items()):
        print(('%9s = %s' % item).replace('\n', '\n' + ' ' * 12))
