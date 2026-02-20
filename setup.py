#!/usr/bin/env python3

import os

import setuptools
from setuptools.command.build_py import build_py as _build_py

class build_py(_build_py):
    """A modified `build_py` command to download data files.
    """

    def run(self):
        # build the rest of the package as normal
        _build_py.run(self)

        # get the path where to download the files
        flexo_blast_path = os.path.join(self.build_lib, "flexo_blast")

setuptools.setup(cmdclass={"build_py": build_py})