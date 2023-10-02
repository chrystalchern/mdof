# Configuration file for the Sphinx documentation builder.
#
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
from pathlib import Path
project = 'mdof'
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
    'css/home-css/'+str(file.name) for file in (Path(__file__).parents[0]/"_static/css/home-css/").glob("vars*.css")
] # + [
#     'css/css/'+str(file.name) for file in (Path(__file__).parents[0]/"_static/css/css/").glob("*.css")
# ]
html_additional_pages = {'index': 'home.html'}
html_context = {
    'description': description,
    'examples': [
        {"title": "SISO Basics",         "link": "examples/01_SISO_Intro", "image": "sdof_full.png"},
        {"title": "SISO Event",          "link": "examples/02_SISO_Event", "image": "examples_02_SISO_Event_20_0.png"},
        {"title": "SIMO Event",          "link": "examples/04_SIMO_Event", "image": "lumps.svg"},
        {"title": "MIMO Basics",         "link": "examples/05_MIMO_Intro", "image": "2dof_full.png"},
        {"title": "MIMO Event",          "link": "examples/05_MIMO_Event"},
        {"title": "SISO History",        "link": "examples/03_SISO_History", "image": "examples_03_SISO_History_7_2.png"},
        {"title": "MIMO History",        "link": "examples/06_MIMO_History"},
        {"title": "ResponseSpectrum",    "link": "examples/ResponseSpectrum"},
        {"title": "transfer function",   "link": "examples/transfer_function"},
        {"title": "juang experiment",    "link": "examples/juang_experiment"},
    ],
    **globals()
}
html_show_sourcelink = False
html_theme_options = {
    "github_url": f"https://github.com/BRACE2/{project}"
}

autodoc_member_order = 'bysource'
