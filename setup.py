#!/usr/bin/env python

"""
Call `pip install -e .` to install package locally for testing.
"""

from setuptools import setup, find_packages


setup(
    name="m6amap",
    version="0.1.0",
    packages=find_packages(),
    #install_requires=[
    #    "pandas",
    #    "numpy",
    #    "gffutils",
    #],
    entry_points={
        "console_scripts": [
            "m6amap=m6amap.__main__:main"
        ]
    },
    author="Selina Chen",
    author_email="yc4635@columbia.edu",
    description="A package to find predicted RNA modification sites (m6A), their associated diseases, genes, and pathways from nanopore sequencing results.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/nya0o0/project",
    classifiers=[
        "Programming Language :: Python :: 3",
        #"License :: OSI Approved :: MIT License",
        #"Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
