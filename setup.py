from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="agrfilegenerator",
    version="3.1.0",
    author="Alliance of Genome Resources",
    author_email="valearna@caltech.edu",
    description="This tool creates files from Alliance resources",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alliance-genome/agr_file_generator",
    include_package_data=True,
    package_dir={'': 'src'},
    packages=find_packages('src'),
    py_modules=['common', 'data_source'],
    install_requires=[
        'neo4j==1.7.3',
        'neobolt==1.7.13',
        'neotime==1.7.4',
        'python-dateutil==2.7.3',
        'requests==2.21.0',
        'Bio==0.1.0',
        'wget==3.2',
        'biopython==1.73',
        'pyfaidx==0.5.5.2',
        'wheel==0.33.4',
        'pytest==4.6.2',
        'Click==7.0',
        'requests==2.21.0',
        'requests-toolbelt==0.9.1',
        'retry==0.9.2',
        'urllib3==1.24.2',
        'coloredlogs==10.0',
        'PyYAML==5.1',
        'json5==0.8.5'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
