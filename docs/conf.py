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
    **globals()
}
html_show_sourcelink = False
html_theme_options = {
    "github_url": f"https://github.com/BRACE2/{project}"
}

autodoc_member_order = 'bysource'
