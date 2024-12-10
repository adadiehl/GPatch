import sys

if sys.version_info > (3, ) and sys.version_info < (3, 7):
    sys.exit("ERROR: GPatch requires Python 3.7 or greater")

try:
    from setuptools import setup, find_packages
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()

from setuptools import setup, find_packages
from pathlib import Path
from glob import glob

def readme():
    if sys.version_info > (3, ):
        with open(Path(__file__).parent.resolve() / 'README.md', encoding='utf-8') as md:
            return md.read()
    else:
        with open('README.md') as md:
            return md.read()

def main():

    metadata = dict(
        name="GPatch",
        version="0.1.0",
        author="Adam Diehl",
        author_email="adadiehl@umich.edu",
        description="Assemble contigs into a chromosome-scalse pseudo-assembly using alignments to a reference sequence.",
        long_description=readme(),
        long_description_content_type="text/markdown",
        url="https://github.com/adadiehl/GPatch",
        packages = find_packages(),
        package_data = {
            'GPatch': [
                "LICENSE",
                "CODE_OF_CONDUCT.md"
            ]
        },
        setup_requires=[
            'Bio',
            'pysam',
            'minimap2'
        ],
        install_requires=[
            'Bio',
            'pysam',
            'minimap2'
        ],
        entry_points={
            'console_scripts': [
                'GPatch.py = GPatch.GPatch:main'
            ]
        },
        classifiers=[
            "Development Status :: 4 - Beta",
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Natural Language :: English"
        ],
        keywords = "genomics, genome assembly",
        include_package_data=True,
        zip_safe=False,
    )
            
    setup(**metadata)
        
if __name__ == "__main__":
    main()
