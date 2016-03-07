#!/usr/bin/env python

from numpy import fromfile, int32, float64
from os.path import abspath, dirname, join
from subprocess import call

def read():
    def MUEB(file, data):
        data['mu*EB'] = fromfile(file, float64, 1)

    def TCMD(file, data):
        data['Tc'] = fromfile(file, float64, 1)

    def TCEB(file, data):
        data['TcEB'] = fromfile(file, float64, 1)

    def EDGE(file, data):
        data['status(Delta0)'] = fromfile(file, int32, 1)
        data['Delta0'] = fromfile(file, float64, 1)

    def IMAG(file, data):
        data['status(Z, Delta)'], n = fromfile(file, int32, 2)

        data['imag. axis'] = {}

        for key in 'omega', 'Z', 'Delta':
            data['imag. axis'][key] = fromfile(file, float64, n)

        data['imag. axis']['phiC'] = fromfile(file, float64, 1)

    def REAL(file, data):
        n = fromfile(file, int32, 1)

        data['real axis'] = {}

        for key in 'omega', 'Re[Z]', 'Im[Z]', 'Re[Delta]', 'Im[Delta]':
            data['real axis'][key] = fromfile(file, float64, n)

    return locals()

read = read()

def load(filename):
    data = {}

    with open(filename, 'rb') as file:
        while True:
            name = file.read(4)

            if name:
                read[name](file, data)
            else:
                break

    return data

def run(executable=join(dirname(abspath(__file__)), 'eb'), **parameters):

    with open('~temporary.in', 'w') as file:

        for parameter, default in [
            ('T', 10.0),
            ('small', 1e-10),
            ('error', 1e-10),
            ('omegaE', 0.02),
            ('lambda', 1.748),
            ('muStar', 0.1),
            ('upper', 10.0),
            ('lower', 5.0),
            ('limit', 100000),
            ('measurable', True),
            ('resolution', 300),
            ('form', 'data'),
            ('epsilon', 1e-15)]:

            print >> file, parameters.get(parameter, default)

    call([executable, '~temporary.in'])

    return load('~temporary.dat')

if __name__ == '__main__':
    print run()
