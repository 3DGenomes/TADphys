from setuptools import setup

setup(
        name         = 'TADdyn',
        version      = '0.1',
        author       = 'Marco Di Stefano, David Castillo',
        author_email = 'marco.di.distefano.1985@gmail.com',
        packages     = ['taddyn', 'taddyn.utils',
                        'taddyn.modelling'],
        platforms = "OS Independent",
        license = "GPLv3",
        description  = 'TADdyn is a Python library that allows to model and explore single or time-series 3C-based data.',
        long_description = (open("README.rst").read()),
        #url          = 'https://github.com/3DGenomes/tadbit',
        #download_url = 'https://github.com/3DGenomes/tadbit/tarball/master',
    )