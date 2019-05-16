from fsac import __version__
from setuptools import setup, find_packages

setup(

    name='fsac',
    version=__version__,
    packages=find_packages(),
    install_requires=[
        'pandas>=0.22.0'
    ],

    author='Dillon Barker',
    author_email='dillon.barker@canada.ca',

    entry_points={
        'console_scripts': [
            'fsac=fsac.main:main'
        ]
    }
)
