from setuptools import setup, find_packages

setup(
    name='ising_model',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'scikit-learn',
        'tqdm'
    ],
)