# CD1bGMM

Reproduces plots for our paper: DeWitt et al., *An MHC-independent T-cell repertoire is enriched during active tuberculosis*, [eLife 2017 🤞]().

## DEPENDENCIES
* scons
* Python 2.7+/3.5+, with modules:
  * scipy
  * matplotlib
  * seaborn
  * pandas
  * biopython
  * scikit-learn
* TCRdist is a submodule, so clone recursively with `git clone --recursive <this repo>`, or if you've already cloned CD1bGMM run `git submodule update --init --recursive` to add the TCRdist submodule.  
    
For installing scons, and python dependencies, [conda](https://conda.io/docs/) is recommended:
```bash
conda install scons scipy matplotlib seaborn pandas biopython scikit-learn
```
A linux spec file for a python2 environment is included (`spec-file.txt`), which may be used to create a conda environment containing all necessary dependencies.
To create an environment called `myenv`, execute `conda create --name myenv --file spec-file.txt`, then activate the environment with `source activate myenv`.

## Generate results

Simply issue the command `scons`, and plot output will be created in a directory named `output`.
