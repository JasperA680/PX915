# Configuration file for the Sphinx documentation builder.

import os
import sys

# Let Sphinx import Python modules from src/python (bare names) and
# from src/ (for modules that use the python.* package prefix internally).
sys.path.insert(0, os.path.abspath("../src/python"))
sys.path.insert(0, os.path.abspath("../src"))

# numpy's character_backward_compatibility_hook asserts that any '=' key in the
# crackfortran AST lives under a 'vars' parent, but character(len=N) length
# selectors also use '=' and sit one level deeper, triggering an AssertionError
# with an empty message.  Remove the hook before sphinxfortran calls crackfortran.
try:
    import numpy.f2py.crackfortran as _crack
    _crack.post_processing_hooks = [
        h for h in _crack.post_processing_hooks
        if h.__name__ != 'character_backward_compatibility_hook'
    ]
except Exception:
    pass

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
    "../src/fortran/ns_model.f90",
    "../src/fortran/junction.f90",
    "../src/fortran/network_builder.f90",
    "../src/fortran/lane_change.f90",
    "../src/fortran/road_network.f90",
    "../src/fortran/vehicle.f90",
    "../src/fortran/network_init.f90",
    "../src/fortran/network_simulation.f90",
    "../src/fortran/network_io.f90",
    "../src/fortran/nc_config.f90",
    "../src/fortran/run_network.f90",
    "../src/fortran/fundamental_diagram.f90",
]

# Use NumPy-style Python docstrings
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']