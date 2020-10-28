__docs__ = """
pip install .
"""

import os, setuptools

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, "README.md"), "r") as fh:
    long_description = fh.read()

requires = [
            'pandas',
            'pyrosetta',
            'beautifulsoup4',
            'lxml'
            ]


setuptools.setup(
    name="pyrosetta_documentarian",
    version="0.0.1",
    author="Matteo Ferla",
    author_email="matteo.ferla@gmail.com",
    description="helper to figure out how to use a pyrosetta mover",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requires,
    url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)