from setuptools import setup, find_packages

setup(
    name='vpk_model',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
    ],
    author='Valery Kubarovsky',
    author_email='vpk@jlab.org',
    description='A package for meson electroproduction calculations and plotting',
    url='https://github.com/yourusername/vpk_model',  # Replace with your actual URL
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
