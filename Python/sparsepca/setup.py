import setuptools

with open('README.md', 'r') as f:
    long_description = f.readlines()
    long_description = ''.join(long_description[4:])

setuptools.setup(
    name="sparsepca",
    version="0.1.1",
    author="Justin Huang and Benjamin Jochem",
    author_email="justinhuang496@gmail.com",
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
