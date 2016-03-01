#!/usr/bin/env python

from numpy import fromfile, int32, float64

def load(filename):
    with open(filename) as file:
        data = {}

        data['status(Z, Delta)'], n = fromfile(file, int32, 2)

        data['imag. axis'] = {}

        for key in 'omega', 'Z', 'Delta':
            data['imag. axis'][key] = fromfile(file, float64, n)

        data['phiC'], data['Tc'] = fromfile(file, float64, 2)

        if file.read(1) == 'T':
            N = fromfile(file, int32, 1)

            data['real axis'] = {}

            for key in 'omega', 'Re[Delta]', 'Im[Delta]':
                data['real axis'][key] = fromfile(file, float64, N)

            data['Delta0'] = fromfile(file, float64, 1)

            data['status(Delta0)'] = fromfile(file, int32, 1)

    return data

if __name__ == '__main__':
    print load('example.dat')
