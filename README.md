# Recovery of 21 cm intensity maps with sparse component separation

> Author: Isabella P. Carucci <a itemprop="sameAs" content="https://orcid.org/0000-0001-5287-0065" href="https://orcid.org/0000-0001-5287-0065" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> <br/> 
> Year: 2020 <br/>
> Email: [ipcarucci@gmail.com](mailto:ipcarucci@gmail.com)
---

**DOI**: [10.1093/mnras/staa2854](http://dx.doi.org/10.1093/mnras/staa2854)   
**arXiv**: [2006.05996](https://arxiv.org/abs/2006.05996)

This is the code repository for <a href="http://dx.doi.org/10.1093/mnras/staa2854" target_="blanck">Recovery of 21 cm intensity maps with sparse component separation</a>. 
All main routines are implemented in Python (>=3.5) and we provide two <a href="https://jupyter-notebook.readthedocs.io/en/stable/" target_="blanck">Jupyter notebooks</a> with demonstrations for reproducing the results presented in the paper.


## Contents

1. [Data](#Data)
1. [Requirements](#Requirements)
1. [Notebooks](#Notebooks)
1. [Credit](#Credit)
1. [License](#License)

### Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3991818.svg)](https://doi.org/10.5281/zenodo.3991818)

The simulated data is provided <a href="https://doi.org/10.5281/zenodo.3991818" target_="blanck">here</a>. 

### Requirements

In order to run the notebooks, the following have to be installed:

* <a href="https://www.python.org/" target_="blank">Python</a> (require >=3.5)
* <a href="http://jupyter.org/" target_="blank">Jupyter</a> (recommend >=1.0.0)
* <a href="https://matplotlib.org/" target_="blank">Matplotlib</a> (recommend >=3.1.1)
* <a href="http://www.numpy.org/" target_="blank">NumPy</a> (recommend >=1.16.4)
* <a href="https://www.scipy.org/" target_="blank">SciPy</a> (recommend >=1.3.0)
* <a href="https://github.com/healpy/" target_="blank">healpy</a> (recommend >=1.12.9)
* <a href="https://www.h5py.org/" target_="blank">h5py</a> (recommend >=2.10.0)


### Notebooks

1. [Prepare the dataset](./notebooks/tut_1_preparation.ipynb)

The objective of this notebook is to provide a first glimpse of the simulation and prepare the data-set for the foreground-removal exercise. After having loaded the simulation, we add white instrumental noise for each channel. We smooth all maps with a frequency-dependent Gaussian beam. We re-smoothed back to the larger beam so all maps share a common resolution.

2. [Perform the BSS](./notebooks/tut_2_perform_GMCA.ipynb)

In this notebook we perform the Blind Source Separation of the maps. We first project the data in wavelet-space and than we apply the [GMCA](http://md.cosmostat.org/Generalized_MCA.html) algorithm, Generalised Morphological Component Analysis [(Bobin et al. 2007)](http://dx.doi.org/10.1109/TIP.2007.906256). We check results by plotting the maps and their power spectra.


### Credit

All credit for original GMCA code goes to [Jérôme Bobin](http://jerome-bobin.fr/).

If you use any of this material in your research, we kindly ask that you cite the <a href="http://dx.doi.org/10.1093/mnras/staa2854" target_="blanck">main paper</a>. 

Regarding bugs, suggestions, questions or feedbacks, feel free to contact the author ([ipcarucci@gmail.com](mailto:ipcarucci@gmail.com)).


### License

GNU General Public License v3.0, read the [LICENSE](LICENSE) for more information.


