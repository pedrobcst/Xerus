


XERUS (X-Ray Estimation and Refinement Using Similarity)
========================================================
![img](img/g163.png)

Paper
====
This is the main repository for
  * Baptista de Castro, P., Terashima, K., Esparza Echevarria, M.G., Takeya, H. and Takano, Y. (2022), XERUS: An Open-Source Tool for Quick XRD Phase Identification and Refinement Automation. Adv. Theory Simul. 2100588. https://doi.org/10.1002/adts.202100588

Introduction
============
Welcome to the _Xerus_ project. Xerus is an open-source python wrapper / plugin around the GSASII Scriptable package,
for automazation of phase quantification and Rietveld analysis by combinining similarity calculation of simulated patterns
(Pearson`s) with quick Rietveld refinements.

_Xerus_ is only possible due to the existence of the following projects
* COD (Crystollographic Open Database)
* The Materials Project (MP)
* AFLOW Database
* OQMD (Open Quantum Materials Database)  
* GSASII Scriptable Engine
* pymatgen

_Xerus_ is designed to perform analysis through __Jupyter__ notebooks, providing an easy to use API from phase quantification to Rietveld optimization

The main mechanisms behind _Xerus_ stems from the following papers:

* __Clustering of XRD patterns using similarity__:
    * Iwasaki, Y., Kusne, A.G. & Takeuchi, I. Comparison of dissimilarity measures for cluster analysis of X-ray diffraction data from combinatorial libraries. npj Comput Mater 3, 4 (2017). https://doi.org/10.1038/s41524-017-0006-2

* __Optmization of Rietveld refinements__:
    * Ozaki, Y., Suzuki, Y., Hawai, T. et al. Automated crystal structure analysis based on blackbox optimisation. npj Comput Mater 6, 75 (2020). https://doi.org/10.1038/s41524-020-0330-9

Plus our own pattern removal iterative process coupled with quick rietveld refinements that allows for multiphase characterization.


Installation
============

In this section we will briefly introduce how to install _Xerus_ in the easiest manner possible.
## Requirements
### Before installing Xerus OS and Services Requeriment:
#### System OS
* _Xerus_ only supports Linux based and macOS systems.
_Xerus_ was mostly developed in __Ubuntu 20.04__ and have been also tested in __CentOS 7.x__ and __MacOS systems (Intel)__ (M1 machines not available for testing.)
```diff
- We have not tested in Windows OS yet.
```

#### Materials Project APIKey
Xerus relies on the Materials Project API for downloading crystal structures with requested chemical space.
Therefore, please register and obtain an suitable API key (Free) at: \
www.materialsproject.org \
After registering, you can check your API Key by clicking the API tab at the upper right side of the website.

##### MongoDB
_Xerus_ relies on a MongoDB server for caching crystal structures downloaded from the providing databases (COD, AFLOW and MP) \
To install the community version (Free) and run please follow the steps listed in:
 * [Ubuntu Installation Steps](https://docs.mongodb.com/manual/tutorial/install-mongodb-on-ubuntu/)
 * [MacOS Installation Steps](https://docs.mongodb.com/manual/tutorial/install-mongodb-on-os-x/)

## _Xerus_
### Package Installation
_Xerus_ currently can only be installed via conda. Therefore, if you do not have conda, please follow the install instructions [here](https://www.anaconda.com/products/individual). \
To install _Xerus_ using a virtual enviroment with Python 3.8 follow this steps:

```
conda create -n Xerus python==3.8 anaconda
conda activate Xerus
git clone http://www.github.com/pedrobcst/Xerus/
cd Xerus
pip install -r requirements.txt
```
> :warning: You might have trouble installing pymatgen if gcc is not present in your system. You can them for example do sudo apt install g++ to install in Ubuntu, then run pip install -r requirements.txt again. 


### Configuration

If all the packages installation are successful, it is needed to correctly set the configuration file at `Xerus/settings/config.conf`:

#### Necessary
* __[mongodb]__
  * _host_: defaults to _localhost_
    * Enter ip address to mongodb server.
  * _user_: username (if necessary)
    * If configuraed to use authentication, please provide the username here
  * _password_: username password (if necessary)
    
* __[mp]__
  * _apikey_: Materials project API key.
    * Please provide your materials project API key here
  
#### Optional
* _Extra_ config:
  * __[gsasii]__ 
    * _instr_params_: Path to .instrprm (GSASII)
      * _Xerus_ by default comes with an .instprm file obtained from fitting NIST Si in our XRD Machine (Rigaku MiniFlex 600).
  It is _recommended_ that you follow the GSASII tutorial for obtaining an .instrprm for your own machine but probably not necessary.
      
    * _testxrd_: 
      * _Xerus_ by default use one of the _Examples_ data for testing purposed of obtained CIFs from the repositories. Feel free to change this to any of your liking

### Testing installation

If all the above steps were done sucessfuly (pip install, mongo running [locally] and API key set),
please do the following:
```bash
cd tests
pytest
```
If all tests sucessfuly pass, _Xerus_ should be ready for use.

> :warning: This process might take a while. 
# Usage

- To learn how to use Xerus, please follow the Notebook, located at Examples/Examples.ipynb

> :warning: Make sure that before running the examples, you have started the MongoDB server. 

# Citing
If you use _Xerus_ please __also__ cite the following papers:

* If you use only phase matching:
  * Baptista de Castro, P., Terashima, K., Esparza Echevarria, M.G., Takeya, H. and Takano, Y. (2022), XERUS: An Open-Source Tool for Quick XRD Phase Identification and Refinement Automation. Adv. Theory Simul. 2100588. https://doi.org/10.1002/adts.202100588
  * Toby, B. H., & Von Dreele, R. B. (2013). "GSAS-II: the genesis of a modern open-source all purpose crystallography software package". Journal of Applied Crystallography, 46(2), 544-549. ​doi:10.1107/S0021889813003531


* If you use blackbox method for refinement please also cite:
  * Ozaki, Y., Suzuki, Y., Hawai, T. et al. Automated crystal structure analysis based on blackbox optimisation. npj Comput Mater 6, 75 (2020). https://doi.org/10.1038/s41524-020-0330-9

# License
The code is licensed under an MIT license.
The data used for benchmarking is licensed under CC4.0 and was taken from

Szymanski, N. J., Bartel, C. J., Zeng, Y., Tu, Q., & Ceder, G. (2021). Probabilistic Deep Learning Approach to Automate the Interpretation of Multi-phase Diffraction Spectra. Chemistry of Materials.
