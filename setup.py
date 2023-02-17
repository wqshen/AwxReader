# -*- coding: utf-8 -*-
# @Author: wqshen
# @Date: 2019/4/6 11:59
# @Last Modified by:   wqshen


from os.path import exists
from setuptools import setup, find_packages

name = 'awx'
version = "0.1"

install_requires = [
    'numpy',
    'xarray',
    'pyproj'
]

classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research/Operation',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Topic :: GeoScience :: AWX Data Reader',
]

setup(
    name=name,
    version=version,
    description='AWX Satellite Data Reader',
    author='wqshen',
    author_email='wqshen91@163.com',
    long_description=open('README.rst').read() if exists('README.rst') else '',
    python_requires='>=3.7',
    install_requires=install_requires,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'awx_info=awx:_info',
            'awx_to_nc=awx:_convert_to_nc',
        ]
    },
    include_package_data=True,
    zip_safe=False,
)
