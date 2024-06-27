import os
import setuptools

setuptools.setup(
    name = "simple_bez",
    version = "0.10",
    author = "Brumo Maximilian Voss",
    author_email = "bruno.m.voss@gmail.com",
    description = ("simple bezier and not quite spline functions"),
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    #license = "MIT",
)
