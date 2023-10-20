# Configuration file for the Sphinx documentation builder.
#
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
from pathlib import Path
project = 'mdof'
copyright = '2023, Chrystal Chern'
author = 'Chrystal Chern'
description = "Fast and friendly structural system identification."
version = '0.0.6'
release = '0.0.6'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    # 'autoapi.extension',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
    'nbsphinx'
]


# autoapi_dirs = ['../src/mdof']


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


source_suffix = '.rst'
root_doc = 'index'
language = 'en'

# -- Options for HTML output -------------------------------------------------
html_title = project
html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']
html_favicon = './_static/images/favicon.ico'
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
        {"title": "Overview",            "link": "examples/00_Overview",   "image": "../_static/images/gallery/overview.svg"},
        {"title": "SISO Basics",         "link": "examples/01_SISO_Intro", "image": "../_static/images/gallery/sdof_full.svg"},
        {"title": "MIMO Basics",         "link": "examples/04_MIMO_Intro", "image": "../_static/images/gallery/2dof_full.svg"},
    ],
    **globals()
}
html_show_sourcelink = False
html_theme_options = {
    "github_url": f"https://github.com/BRACE2/{project}"
}

autodoc_member_order = 'bysource'
