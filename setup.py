from distutils.core import setup
setup(name='pv_atmos',
      version='2.3.2',
      description='Utilities for scientific visualization with ParaView',
      long_description='Adds convenient functionality to visualize netCDF data in arbitrary coordinates, including logarithmic and spherical geometries. Mainly created for atmospheric and oceanic sciences in mind, it is general enough to be useful for many other applications. Read more on the GitHub page, or in the peer-reviewed open access article, which can be found at http://dx.doi.org/10.5334/jors.al',
      author='Martin Jucker',
      author_email='coding@martinjucker.com',
      license='MIT',
      url='https://github.com/mjucker/pv_atmos',
      package_dir={'pv_atmos': ''},
      packages=['pv_atmos'],
      )
