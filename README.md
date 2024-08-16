# DeepEC
This is a fork of [DeepEC](https://bitbucket.org/kaistsystemsbiology/deepec/src/master/). Changes were made for integration with newer versions of [E2P2](https://github.com/carnegie/E2P2).

## Procedure

**Note**: 
Size of the protein sequence input file should be adjusted according to the memory size of your computer. 
This source code was developed in Linux, and has been tested in Ubuntu 14.04 with Python 2.7, Python 3.4, Python 3.5 or Python 3.6. 
It should be noted that Python 3.7 is currently not supported.

1. Clone the repository

        git clone git@github.com:bxuecarnegie/deepec.git

2. Create and activate a conda environment

        conda env create -f environment.yml
        conda activate deepec

## Example

- Run DeepEC

        python deepec.py -i ./example/test.fa -o ./output 

