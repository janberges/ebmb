#!/usr/bin/env python

from numpy import fromfile, int32, float64
from os.path import abspath, dirname, join
from subprocess import call

def load(filename):
    with open(filename, 'rb') as file:
        data = {}

        data['status(Z, Delta)'], n = fromfile(file, int32, 2)

        data['imag. axis'] = {}

        for key in 'omega', 'Z', 'Delta':
            data['imag. axis'][key] = fromfile(file, float64, n)

        data['phiC'], data['mu*EB'], data['Tc'] = fromfile(file, float64, 3)

        if file.read(1) == 'T':
            N = fromfile(file, int32, 1)

            data['real axis'] = {}

            for key in 'omega', 'Re[Delta]', 'Im[Delta]':
                data['real axis'][key] = fromfile(file, float64, N)

            data['Delta0'] = fromfile(file, float64, 1)

            data['status(Delta0)'] = fromfile(file, int32, 1)

    return data

def run(executable=join(dirname(abspath(__file__)), 'eb'), **parameters):

    with open('~temporary.in', 'w') as file:

        for parameter, default in [
            ('T', 10.0),
            ('omegaE', 0.02),
            ('lambda', 1.748),
            ('muStar', 0.1),
            ('upper', 0.2),
            ('lower', 0.1),
            ('continue', True),
            ('resolution', 300),
            ('limit', 100000),
            ('tiny', 1e-15),
            ('form', 'data')]:

            print >> file, parameters.get(parameter, default)

    call([executable, '~temporary.in'])

    return load('~temporary.dat')

if __name__ == '__main__':
    print run()
