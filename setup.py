from setuptools import setup
from distutils.core import setup, Extension

def main():
    # c++ module to compute the distance matrix of single model
    squared_distance_matrix_module = Extension('taddyn.squared_distance_matrix',
                                               language = "c++",
                                               runtime_library_dirs=['3d-lib/'],
                                               sources=['src/3d-lib/squared_distance_matrix_calculation_py.c'],
                                               extra_compile_args=["-ffast-math"])

    setup(
        name         = 'TADdyn',
        version      = '0.1',
        author       = 'Marco Di Stefano, David Castillo',
        author_email = 'marco.di.distefano.1985@gmail.com',
        ext_modules  = [squared_distance_matrix_module],
        packages     = ['taddyn', 'taddyn.utils',
                        'taddyn.modelling'],
        platforms = "OS Independent",
        license = "GPLv3",
        description  = 'TADdyn is a Python library that allows to model and explore single or time-series 3C-based data.',
        long_description = (open("README.rst").read()),
        #url          = 'https://github.com/3DGenomes/tadbit',
        #download_url = 'https://github.com/3DGenomes/tadbit/tarball/master',
    )

if __name__ == '__main__':

    exit(main())
