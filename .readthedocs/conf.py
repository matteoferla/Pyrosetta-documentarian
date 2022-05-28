# sphinx-apidoc -o .readthedocs . sphinx-apidoc --full -A 'Matteo Ferla';

project = 'pyrosetta-help'
author = 'Matteo Ferla'
copyright = '2022, University of Oxford'
github_username = 'matteoferla'
github_repository = 'pyrosetta_help'

extensions = [
    'readthedocs_ext.readthedocs',
    #'sphinx_toolbox.more_autodoc'  test required
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.todo',
]

templates_path = ['_templates']

language = 'en'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


def skip(app, what, name, obj, would_skip, options):
    if name in ( '__init__',):
        return False
    return would_skip

def setup(app):
    app.connect('autodoc-skip-member', skip)

html_theme = 'sphinx_rtd_theme'
html_static_path = []  # '_static'

todo_include_todos = True

# =============================================================================

import m2r2  # noqa
import os, re
from typing import Optional

repo_base_path = os.path.abspath("../")

def convert_write(markdown_filename, srt_filename, change_title:Optional[str]=None):
    # unlike Fragmenstein there are no images to convert
    # so we can just copy the file
    with open(markdown_filename) as fh:
        markdown_block = fh.read()
    if change_title:
        markdown_block = re.sub(r'#+\s*(.*)', r'# ' + change_title, markdown_block, 1)
    #markdown_block = re.sub(r'\[(?P<label>.*?)\]\((?P<link>.*?)\)', fix_md_link, markdown_block)
    rst_block = m2r2.convert(markdown_block)
    with open(srt_filename, 'w') as fh:
        fh.write(rst_block)

convert_write(os.path.join(repo_base_path, 'README.md'), 'introduction.rst', 'Overview')
