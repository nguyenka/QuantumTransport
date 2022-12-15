#!/usr/bin/env python3
# pip install --user -e .
from setuptools import setup   

requirements = ['ase>=3.22.1','cclib>=1.7.1', 'matplotlib', 'munch>=2.5.0','numpy>=1.19.5']

setup(name="QuantumTransport",
      version="0.1",
      description="Calculating electron spin polarization and transmission",
      author="Kiet A. Nguyen",
      author_email="kiet.nguyen.2.ctr@afrl.af.mil",
      url="https://github.com/nguyenka/QuantumTransport",
      license="GNU General Public License 3.0",
      platforms="platform independent",
      packages=[ "QuantumTransport", ],
      scripts=[
          "bin/qt",
          "bin/qt_dat",
          "bin/g4qt",
          "bin/t4qt"
          ],
      install_requires=requirements,
      classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.8",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        ],
      long_description="""
Calculating electron spin polarization and transmission
-------------------------------------------------------

This package calculates spin polarization and transmission for:
1) Chiral induced spin polarization transmission with SOC a closed-shell singlet.
1) Spin polarization transmission for an open shell radical.
""")

