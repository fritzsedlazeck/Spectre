from setuptools import setup, find_packages

setup(
    name='spectre-cnv',
    version='0.2.0',
    packages=find_packages(),
    url='https://github.com/fritzsedlazeck/Spectre',
    license='MIT',
    author='Philippe Sanio',
    author_email='philippe.sanio@gmail.com',
    description='Long read copy number variation (CNV) caller',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
