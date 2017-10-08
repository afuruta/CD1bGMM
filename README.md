# CD1bGMM

Reproduces plots for our paper: DeWitt et al., *An MHC-independent T-cell repertoire is enriched during active tuberculosis*, [eLife 2017 🤞]().

## DEPENDENCIES
* scons 3+
* Python 3.5+, with modules:
  * scipy
  * matplotlib
  * seaborn
  * pandas
  * biopython
  * scikit-learn

**NOTE:** for installing scons, and python dependencies, [conda](https://conda.io/docs/) is recommended:
```bash
conda install scons scipy matplotlib seaborn pandas biopython scikit-learn
```
Alternatively, an example linux environment spec file is included (`spec-file.txt`), which may be used to create a conda environment.
For example, to create an environment called `myenv`, execute `conda create --name myenv --file spec-file.txt`, then activate the environment with `source activate myenv`.
