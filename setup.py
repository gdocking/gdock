from setuptools import setup, find_packages

setup(
    name="gdock",
    license="BSD Zero Clause (0BSD)",
    version="1.1.6",
    author="Rodrigo V. Honorato",
    description="Genetic Algorithm applied to Protein-Protein Docking",
    author_email="rvhonorato@protonmail.com",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Zero-Clause BSD (0BSD)",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.9, <4",
    install_requires=[
        "toml",
        "biopython",
        "scipy",
        "pandas",
        "deap==1.3.1",
        "mgzip",
        "pdb-tools>=2.4.1",
    ],
    entry_points={
        "console_scripts": [
            "gdock=gdock.cli:main",
        ],
    },
)
