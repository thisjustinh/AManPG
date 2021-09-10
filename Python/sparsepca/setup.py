import setuptools

with open('README.md', 'r') as f:
    readme = f.readlines()[4:]
    readme = ''.join(readme)

with open('HISTORY.md', 'r') as f:
    history = f.read()

setuptools.setup(
    name="sparsepca",
    version="0.2.3",
    author="Shixiang Chen, Justin Huang, Benjamin Jochem, Shiqian Ma, Lingzhou Xue, and Hui Zou",
    author_email="lzxue@psu.edu",
    description="Sparse Principal Component Analysis in Python",
    long_description=''.join([readme, '\n\n', history]),
    long_description_content_type="text/markdown",
    # packages=setuptools.find_packages(),
    license='MIT',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    project_urls={
        'Documentation' : 'https://xinging-birds.github.io/AManPG/',
        'Source': 'https://github.com/xinging-birds/AManPG',
    },
    python_requires=">=3.8",
    py_modules=["sparsepca"],
    install_requires=["numpy"]
)
