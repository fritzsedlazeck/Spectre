from setuptools import setup, find_packages

setup(name="spectre",
        version="",
        description="",
        author="",
        url="",
        include_package_data=True,
        packages=find_packages(exclude=[]),
        python_requires='>=3.7',
        install_requires=["pysam>=0.21.0", "numpy>=1.24.3", "pandas>=2.0.1", "matplotlib>=3.7.1", "scipy>=1.10.1"],
        long_description="",
        classifiers=[],
        platforms=["any"]
        )
