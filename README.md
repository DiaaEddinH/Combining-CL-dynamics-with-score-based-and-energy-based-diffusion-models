# Combining complex Langevin dynamics with score-based and energy-based diffusion models

This repository contains the workflow and data used in [Combining complex Langevin dynamics with score-based and energy-based diffusion models.](https://arxiv.org/abs/2510.01328) 


## Requirements
* Conda such as from [miniforge](https://docs.conda.io/en/latest/miniconda.html)
* Python 3.11 or later
* Snakemake


## Setup
1. Set an environment with the above requirements
2. Clone repository and `cd` to it:
```gitbash
git clone --recurse-submodules https://github.com/DiaaEddinH/Combining-CL-dynamics-with-score-based-and-energy-based-diffusion-models.git

```

3. Run pip install

```shellbash
pip install -e libs/DiffusionModels
```

4. Unzip data files `samples.zip` & `weights.zip` from Zenodo and place them in `data` directory.

Alternatively, you may unzip the `Archive.zip` file from and run steps 2-4.

The data and workflow can be found in 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17725665.svg)](https://doi.org/10.5281/zenodo.17725665)

Assets used by the publication are obtained by running snakemake:

```Snakemake
snakemake --cores <number_of_cores> --config EXPERIMENT=QMODEL|GMODEL --use-conda
```
