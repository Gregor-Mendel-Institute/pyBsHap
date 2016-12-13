from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pyBsHap',
    version='0.0.1',
    description='refined bs-seq toolkit',
    long_description=long_description,
    url='https://github.com/rahbz/pyBsHap',
    author=['Rahul Pisupati'],
    author_email='rahul.bharadwaj.p@gmail.com',
    license='GMI',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.4',
    ],
    keywords='bisulfite sequencing',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=[
        "numpy >=1.9.0",
        "pysam",
        "pyfaidx"
    ],
    entry_points={
        'console_scripts': [
            'bshap=bshap:main'
        ],
    },
)
