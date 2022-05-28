__docs__ = """
pip install .
"""

import os, setuptools

here = os.path.abspath(os.path.dirname(__file__))
if os.path.exists("README.md"):
    with open("README.md", "r") as fh:
        long_description = fh.read()
else:
    long_description = 'see [GitHub Repo](https://github.com/matteoferla/Pyrosetta-documentarian)'

if os.path.exists("requirements.txt"):
    with open("requirements.txt", "r") as fh:
        requirements = [line.split("#")[0].strip() for line in fh.readlines()]
else:
    requirements = []

setuptools.setup(
    name="pyrosetta_documentarian",
    version="0.0.1",
    author="Matteo Ferla",
    author_email="matteo.ferla@gmail.com",
    description="helper to figure out how to use a pyrosetta mover",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=requirements,
    include_package_data=True,
    url="https://github.com/matteoferla/Pyrosetta-documentarian",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)