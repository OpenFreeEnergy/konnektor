# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
from pathlib import Path
import shutil
import sys

sys.path.insert(0, os.path.abspath("../src"))


# -- Project information -----------------------------------------------------

project = "konnektor"
copyright = "2022, The OpenFE Development Team"
author = "The OpenFE Development Team"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    "nbsphinx",
    "nbsphinx_link",
    "sphinx.ext.autosectionlabel",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# autodoc_mock_imports = ['lomap', 'networkx', 'openff', 'openff.toolkit', 'rdkit', 'pytest',
#                        'typing_extensions',
#                        'click', 'plugcli']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_static_path = ["_static"]
html_theme = "ofe_sphinx_theme"
autoclass_content = "both"
html_favicon = "_static/img/logo/konnektor_logo_style1.ico"
html_logo = "_static/img/logo/konnektor_logo_style2_borderless.png"
html_theme_options = {
    "logo": {
        "text": "Konnektor Documentation",
    },
    "icon_links": [
        {
            "name": "Github",
            "url": "https://github.com/OpenFreeEnergy/konnektor",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        }
    ],
    "accent_color": "FeelingSick",
}
# html_logo = "_static/Squaredcircle.svg"


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_css_files = [
    "css/custom.css",
]

# Extensions for the myst parser
myst_enable_extensions = [
    "dollarmath",
    "colon_fence",
    "smartquotes",
    "replacements",
    "deflist",
    "attrs_inline",
]
myst_heading_anchors = 3


example_notebooks_path = Path("ExampleNotebooks")

try:
    if example_notebooks_path.exists():
        pass
    else:
        source = Path("../examples")
        shutil.copytree(source, example_notebooks_path)
except Exception as e:
    raise OSError("Could not copy over example notebooks")
