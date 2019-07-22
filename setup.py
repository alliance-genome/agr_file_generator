import os
from setuptools import find_packages, setup

here = os.path.abspath(os.path.dirname(__file__))

def _read_file(name):
    with open(os.path.join(here, name)) as fp:
        return fp.read()

README = _read_file('README.rst')
CHANGES = _read_file('changelog.rst')
INSTALL_REQUIRES = _read_file('requirements.txt').splitlines()

setup(
    name='agr_file_generator',
    version='0.1.dev0',
    url='http://www.alliance.org/',
    author='Adam Wright, Matt Russell',
    author_email='developers@alliancegenome.org',
    description='Variant Call Format file (VCF) genrator for the Alliance Genome Resources',
    long_description="""\
    Provides command line facilities to generate VCF files. 
    """,
    license='MIT',
    keywords='AGR, Model Organisms',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=INSTALL_REQUIRES,
    entry_points={
        'console_scripts': [
            'agr_file_generator=agr.app:main',
        ],
        # TOOD: add hooks for docs generation?
        # 'zest.releaser.releaser.after': [
        #     'publish_azanium_docs_to_github_pages=hooks:deploy_release'
        # ]
    },
    zip_safe=False,
    classifiers=[
        'Development Status :: 1 - Pre-Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Software Development :: Libraries :: Python Modules :: Tools',
    ],
    extras_require={
        'dev': [
            'Sphinx==1.4.3',
            'ghp-import==0.4.1',
            'sphinx_rtd_theme==0.1.9',
            'zest.releaser[recommended]==6.6.4',
            'pytest==4.6.2',
        ]
    }
)
