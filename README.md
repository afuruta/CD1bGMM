# CD1bGMM

This repo reproduces analysis and plots for our paper:

**A diverse lipid antigen-specific T cell receptor repertoire is clonally expanded during active tuberculosis**. DeWitt WS, Yu KKQ, Wilburn DB, Sherwood A, Vignali M, Day CL, Scriba TJ, Robins HS, Swanson WJ, Emerson RO, Bradley PH, Seshadri C. In press at _Journal of Immunology_.

Data: http://dx.doi.org/doi:10.21417/B7QG66

## Dependencies
* scons
* Python 2.7+/3.5+, with modules:
  * scipy
  * matplotlib
  * seaborn
  * pandas
  * biopython
  * scikit-learn
* TCRdist is a submodule, so clone recursively when you clone this repo with  
  `git clone --recursive https://github.com/seshadrilab/CD1bGMM.git`,  
  or if you've already cloned this repo, run  
  `git submodule update --init --recursive` to add the TCRdist submodule.  

For installing scons, and python dependencies, [conda](https://conda.io/docs/) is recommended:
```bash
conda install scons scipy matplotlib seaborn pandas biopython scikit-learn
```
Alternatively, a linux spec file for a python2 environment is included (`spec-file.txt`), which may be used to create a conda environment containing all necessary dependencies.
To create an environment called `myenv`, execute  
`conda create --name myenv --file spec-file.txt`,  
then activate the environment with  
`source activate myenv`.

## Generate results

Simply issue the command `scons`. Necessary data will be downloaded, and plot output will be created in a directory named `output`.
