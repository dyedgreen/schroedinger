# Build the Numerov c module for use in local Python
# (This does NOT install the module, but generates
# a local file that allows the module to be called.)

from distutils.core import setup, Extension

sources = [
  "module.c",
  "eigenvalues.c",
  "numerov.c",
]

module_numerov = Extension("numerov_c", sources)

setup (
  version = '0.1',
  description = 'The c implementation of the Numerov module',
  ext_modules = [module_numerov]
)