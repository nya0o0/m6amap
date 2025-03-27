#!/usr/bin/env python

"""
Call `pip install -e .` to install package locally for testing.
"""

from setuptools import setup, find_packages


setup(
    name="m6alinker",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "gffutils",
    ],
    entry_points={
        "console_scripts": [
            "m6alinker=m6alinker.__main__:main"
        ]
    },
    author="Selina Chen",
    author_email="yc4635@columbia.edu",
    description="A package to process m6A modification sites with GTF annotation, providing details for searching in M6ADD database.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/nya0o0/M6A-Linker",
    classifiers=[
        "Programming Language :: Python :: 3",
        #"License :: OSI Approved :: MIT License",
        #"Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)