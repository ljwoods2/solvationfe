"""
Location of data files
======================

Use as ::

    from solvationfe.data.files import *

"""

__all__ = [
    "MDANALYSIS_LOGO",  # example file of MDAnalysis logo
]

import importlib.resources

data_directory = importlib.resources.files("solvationfe") / "data"

MDANALYSIS_LOGO = data_directory / "mda.txt"
