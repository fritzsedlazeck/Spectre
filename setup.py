from setuptools import setup, find_packages

setup(name="spectre",
        version="0.1.0",
        description="Spectre - a CNV caller for long read data",
        author="",
        url="https://github.com/fritzsedlazeck/Spectre",
        include_package_data=True,
        packages=find_packages(exclude=[]),
        entry_points = {
                'console_scripts': [
                        'spectre = spectre.spectre:main',                  
                ],              
        },
        python_requires='>=3.7',
        install_requires=["pysam>=0.21.0", "numpy>=1.24.3", "pandas>=2.0.1", "matplotlib>=3.7.1", "scipy>=1.10.1"],
        long_description="",
        classifiers=[
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Natural Language :: English",
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10"
        ],
        platforms=["any"]
        )
