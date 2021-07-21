import setuptools

with open('README.md', 'r') as f:
    readme = f.read()

with open('HISTORY.md', 'r') as f:
    history = f.read()

setuptools.setup(
    name="sparsepca",
    version="0.2.1",
    author="Justin Huang, Benjamin Jochem, Shiqian Ma, and Lingzhou Xue",
    author_email="lzxue@psu.edu",
    description="Sparse Principal Component Analysis in Python",
    long_description=''.join([readme, '\n\n', history]),
    long_description_content_type="text/markdown",
    # packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    py_modules=["sparsepca"],
    install_requires=["numpy"]
)
