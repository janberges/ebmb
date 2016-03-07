#!/usr/bin/env python

from numpy import fromfile, int32, float64
from os.path import abspath, dirname, join
from subprocess import call

def read():
    def MUEB(file):
        return fromfile(file, float64, 1)[0]

    TCMD = TCEB = MUEB

    def EDGE(file):
        return {
            'status': fromfile(file, int32, 1),
            'Delta0': fromfile(file, float64, 1)}

    def IMAG(file):
        data = {}

        data['status'], n = fromfile(file, int32, 2)

        for key in 'omega', 'Z', 'Delta':
            data[key] = fromfile(file, float64, n)

        data['phiC'], = fromfile(file, float64, 1)

        return data

    def REAL(file):
        n, = fromfile(file, int32, 1)

        return {key: fromfile(file, float64, n)
            for key in ('omega', 'Re[Z]', 'Im[Z]', 'Re[Delta]', 'Im[Delta]')}

    return locals()

read = read()

def load(filename):
    data = {}

    with open(filename, 'rb') as file:
        while True:
            name = file.read(4)

            if name:
                data[name] = read[name](file)
            else:
                break

    return data

def run(executable=join(dirname(abspath(__file__)), 'eb'), **parameters):

    with open('~temporary.in', 'w') as file:

        for parameter, default in [
            ('T', 10.0), # temperature (K)

            ('small', 1e-10), # negligible gap (eV)
            ('error', 1e-10), # error of critical temperature (K)

            ('omegaE', 0.020), # Einstein frequency (eV)
            ('lambda', 1.748), # electron-phonon coupling
            ('muStar', 0.100), # Coulomb pseudo-potential

            ('upper', 10.0), # overall cutoff frequency (eV)
            ('lower',  5.0), # Coulomb cutoff frequency (eV)

            ('limit', 100000), # maximum number of fixed-point steps

            ('measurable', True), # find measurable gap?
            ('resolution',  300), # real axis resolution

            ('form', 'data'), # output format

            ('epsilon', 1e-15)]: # negligible float difference

            print >> file, parameters.get(parameter, default)

    call([executable, '~temporary.in'])

    return load('~temporary.dat')

if __name__ == '__main__':
    print run()
