# AUTHOR: CJC
# CREATE TIME: 2022/11/17
# DESCRIPTION: SIMPLE SCRIPT


"""
setup.py adapted from https://github.com/kennethreitz/setup.py
"""
from setuptools import setup
import setuptools


setup(
    name='pythf',
    version='1.0',
    description='auxiliary script that working on using TOUGH+HYDRATE',
    author='CJC',
    author_email='2020205219@tju.edu.cn',
    python_requires='3.10',
    url='',
    packages=setuptools.find_packages(),
    scripts=['PYTHF'],
    install_requires=['numpy',
                      'pandas',
                      'matplotlib',
                      'fortranformat',
                      'palettable'],
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Microsoft :: Windows :: Windows 10",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Physics"
    ]
)