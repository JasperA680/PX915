# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PX915 Traffic Flow Modelling'
copyright = '2026, Jasper Allen, Tristan McCarthy, Stephan Gambart, Lucas Belz-Koeling'
author = 'Jasper Allen, Tristan McCarthy, Stephan Gambart, Lucas Belz-Koeling'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions =['sphinxfortran.fortran_domain','sphinxfortran.fortran_autodoc']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

fortran_src = ["../src/fortran/pde_*.f90", "../src/fortran/tasep.f90",
               "../src/fortran/simulation.f90", "../src/fortran/fundamental_diagram.f90"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
