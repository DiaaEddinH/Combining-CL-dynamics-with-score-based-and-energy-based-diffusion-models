# Score-based diffusion models


## Setting up environment

# Requirements
* Conda such as from [miniforge](https://docs.conda.io/en/latest/miniconda.html)
* Python 3.11 or later
* Snakemake


# Setup
1. Set an environment with the above requirements
2. Clone repository and `cd` to it:
```gitbash
git clone https://github.com/DiaaEddinH/Combining-CL-dynamics-with-score-based-and-energy-based-diffusion-models.git

```
3. Run pip install

```shellbash
pip install -e libs/DiffusionModels
```

4. Unzip data files `samples.zip` & `weights.zip` from Zenodo and place them in `data` directory.

5. Run the snakefile:

```Snakemake
snakemake --cores <number_of_cores> --config EXPERIMENT=QMODEL|GMODEL --use-conda
```