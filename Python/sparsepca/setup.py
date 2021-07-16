import setuptools

with open('README.md', 'r') as f:
    long_description = f.readlines()
    long_description = ''.join(long_description[4:])

setuptools.setup(
    name="sparsepca",
    version="0.2.0",
    author="Justin Huang, Benjamin Jochem, Shiqian Ma, and Lingzhou Xue",
    author_email="lzxue@psu.edu",
    description="Sparse Principal Component Analysis in Python",
    long_description=long_description,
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
