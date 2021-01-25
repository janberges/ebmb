#!/usr/bin/env python3

# Copyright (C) 2021 Jan Berges
# This program is free software under the terms of the GNU GPLv3 or later.

import setuptools

with open('README.md', 'r', encoding='utf-8') as README:
    long_description = README.read()

setuptools.setup(
    name                          = 'ebmb',
    version                       = '2021.1',
    author                        = 'Jan Berges',
    author_email                  = '',
    description                   = 'Solve multiband Eliashberg equations',
    long_description              = long_description,
    long_description_content_type = 'text/markdown',
    url                           = 'https://bitbucket.org/berges/ebmb',
    py_modules                    = ['ebmb'],
    python_requires               = '>=2.7',
    install_requires              = [],
    classifiers                   = [
        'Programming Language :: Python',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX :: Linux',
        ],
    )
