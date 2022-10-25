import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="parm",
    version="0.1.0",
    python_requires=">=3.8",
    extras_require=[
        "pysb>=1.13.2",
        "cython>=0.29.25",
    ],
    author="Blake A. Wilson",
    author_email="blake.wilson@utdallas.edu",
    description="PAR2 Activation and calcium signaling Reaction Model .",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NTBEL/PARM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
