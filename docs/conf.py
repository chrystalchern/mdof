# Configuration file for the Sphinx documentation builder.
#
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
from pathlib import Path
project = 'ssid'
copyright = '2023, Chrystal Chern'
author = 'Chrystal Chern'
description = "Fast and friendly structural system identification."
version = '0.0.4'
release = '0.0.4'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    'nbsphinx'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


source_suffix = '.rst'
root_doc = 'index'
language = 'en'

# -- Options for HTML output -------------------------------------------------
html_title = project
html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']
html_favicon = '_static/favicon.ico'
html_css_files = [
    "css/peer.css",
] + [
    'css/css/'+str(file.name) for file in (Path(__file__).parents[0]/"_static/css/css/").glob("*.css")
]
html_additional_pages = {'index': 'home.html'}
html_context = {
    'description': description,
    'examples': [
        {"title": "SISO Intro",          "path": "examples/01_SISO_Intro"},
        {"title": "MIMO Event",          "path": "examples/05_MIMO_Event"},
        {"title": "SISO Event",          "path": "examples/02_SISO_Event"},
        {"title": "MIMO Intro",          "path": "examples/04_MIMO_Intro"},
        {"title": "MIMO Intro",          "path": "examples/05_MIMO_Intro"},
        {"title": "SISO History",        "path": "examples/03_SISO_History"},
        {"title": "SIMO Event",          "path": "examples/04_SIMO_Event"},
        {"title": "MIMO History",        "path": "examples/06_MIMO_History"},
        {"title": "ResponseSpectrum",    "path": "examples/ResponseSpectrum"}, 
        {"title": "transfer function",   "path": "examples/transfer_function"},        
        {"title": "juang experiment",    "path": "examples/juang_experiment"},
    ],
    **globals()
}
html_show_sourcelink = False
html_theme_options = {
    "github_url": f"https://github.com/BRACE2/{project}"
}

autodoc_member_order = 'bysource'
