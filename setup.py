from setuptools import setup,find_packages

setup(
        name='misa',    # This is the name of your PyPI-package.
        version='1.0.0',    # Update the version number for new releases
        scripts=['run_misa.py',], # The name of your scipt, and also the command you'll be using for calling it
        description='MISA: a mixed samples analysis tool',
        long_description='MISA stands for MIxed Sample Analysis tool and addresses the problem of phylogenetic placement of DNA sequences of mixed genomes and hybrids into an\already existing reference tree.',
        long_description_content_type='text/plain',
        url='https://github.com/balabanmetin/misa',
        author='Metin Balaban',
        author_email='balaban@ucsd.edu',
        packages=find_packages(),
        zip_safe = False,
        install_requires=['numpy','treeswift','scipy'],
        include_package_data=True
)
