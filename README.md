# Recovery of 21 cm intensity maps with sparse component separation
---

> Author: [Isabella P. Carucci](http://orcid.org/0000-0001-5287-0065)<br/>
> Year: 2020 <br/>
> Email: [ipcarucci@gmail.com](mailto:ipcarucci@gmail.com)


This is the code repository for <a href="https://arxiv.org/abs/2006.05996" target_="blanck">Recovery of 21 cm intensity maps with sparse component separation</a>. 
All main routines are implemented in Python (>=3.5) and we provide two <a href="https://jupyter-notebook.readthedocs.io/en/stable/" target_="blanck">Jupyter notebooks</a> with demonstrations for reproducing the results presented in the paper.


## Contents

1. [Data](#Data)
1. [Requirements](#Requirements)
1. [Notebooks](#Notebooks)
1. [Credit](#Credit)
1. [License](#License)

### Data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3991818.svg)](https://doi.org/10.5281/zenodo.3991818)

The simulated data is provided here. . . 

### Requirements

In order to run the notebooks, the following have to be installed:

* <a href="https://www.python.org/" target_="blank">Python</a> (require >=3.5)
* <a href="http://jupyter.org/" target_="blank">Jupyter</a> (recommend >=1.0.0)
* <a href="https://matplotlib.org/" target_="blank">Matplotlib</a> (recommend >=3.0.3)
* <a href="http://www.numpy.org/" target_="blank">NumPy</a> (recommend >=1.16.2)
MISSING ....

### Notebooks

1. [Prepare the dataset](./tut_1_preparation.ipynb)

The objective of this notebook is to provide a first glimpse at the simulation and prepare the dataset for the foreground-removal exercise. After having loaded the simulation, we add white instrumental noise for each channel. We smooth all maps with a frequency-dependent Gaussian beam. We re-smoothed back to the larger beam so all maps share a common resolution.

2. [Perform the BSS](./tut_2_perform_GMCA.ipynb)

In this notebook we perform the Blind Source Separation of the maps, making using of the wavelet-transformed dataset and the GMCA algorithm. We check results by plotting the maps and their power spectra.


### Credit

All credit for original GMCA code goes to [Jérôme Bobin](http://jbobin.cosmostat.org).

If you use any of this material in your research, we kindly ask that you cite ...

And for the simulated dataset ... 

### License:

??
