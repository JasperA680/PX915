# Configuration file for the Sphinx documentation builder.

import os
import sys

# Let Sphinx import Python modules from src/python
sys.path.insert(0, os.path.abspath("../src/python"))

# -- Project information -----------------------------------------------------

project = 'PX915 Traffic Flow Modelling'
copyright = '2026, Jasper Allen, Tristan McCarthy, Stephan Gambart, Lucas Belz-Koeling'
author = 'Jasper Allen, Tristan McCarthy, Stephan Gambart, Lucas Belz-Koeling'
release = '0.1.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinxfortran.fortran_domain',
    'sphinxfortran.fortran_autodoc',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

fortran_src = [
    "../src/fortran/pde_*.f90",
    "../src/fortran/tasep.f90",
    "../src/fortran/simulation.f90",
    "../src/fortran/fundamental_diagram.f90",
    "../src/fortran/ns_model.f90",
]

# Use NumPy-style Python docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']