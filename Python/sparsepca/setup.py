import setuptools

with open('README.md', 'r') as f:
    long_description = f.read()

setuptools.setup(
    name="sparsepca",
    version="1.0.0",
    author="Justin Huang and Benjamin Jochem",
    description="Sparse Principal Component Analysis in Python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    py_modules=["sparsepca"],
    package_dir={'':'sparsepca/src'},
    install_requires=[]
)