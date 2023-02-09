import setuptools

# with open("README.md", encoding="utf8") as fh:
#     long_description = fh.read()

setuptools.setup(
    name="parm",
    version="0.3.0",
    python_requires=">=3.8",
    install_requires=["pysb>=1.13.2"],
    extras_require={"cython": "cython>=0.29.25"},
    author="Blake A. Wilson",
    author_email="blake.wilson@utdallas.edu",
    description="PAR2 Activation and calcium signaling Reaction Model.",
    # long_description=long_description,
    # long_description_content_type="text/markdown",
    url="https://github.com/NTBEL/PARM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
