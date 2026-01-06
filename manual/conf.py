# Configuration file for the Sphinx documentation builder.

# https://www.sphinx-doc.org/en/master/usage/configuration.html

project = 'ebmb'
copyright = '2016-2026 Jan Berges'
author = 'Uni Bremen'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'numpydoc',
    'myst_parser',
]

html_theme = 'sphinx_rtd_theme'
html_logo = '../logo/ebmb.svg'
html_theme_options = {
    'logo_only': True,
    'style_nav_header_background': '#e7f2fa',
}
html_css_files = ['style.css']
html_static_path = html_css_files
