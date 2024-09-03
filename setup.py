from setuptools import find_packages, setup

setup(
    name='swarm_interaction_network',
    version='1.0.0',
    description='Here we assess swarm intelligence techniques by examining the social interactions within the swarm.',
    author='Marcos Oliveira',
    url='https://github.com/macoj/swarm_interaction_network.git',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'pandas',
        'networkx',
        'scipy',
        'igraph',
        'fastcluster'
    ],
    include_package_data=True,
    license='GPL3',
)
