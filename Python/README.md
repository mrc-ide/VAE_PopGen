# VAE_PopGen

## Conda environment

We recommend setting up a [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment.

```
conda create --name popgen python=3.12
conda activate popgen
pip install "jax[cpu]"
conda install conda-forge::numpyro
pip install git+ssh://git@github.com/MLGlobalHealth/dsp.git
pip install git+ssh://git@github.com/MLGlobalHealth/sps.git
conda install conda-forge::arviz
conda install hydra-core
conda install jupyter
```

If working on a machine with GPU, the lines installing `jax` and `numpyro` need to be replaced. The full set of commands is as follows: 

```
conda create --name popgen python=3.12
conda activate popgen
pip install --upgrade "jax[cuda12]"
pip install 'numpyro[cuda]' -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
pip install git+ssh://git@github.com/MLGlobalHealth/dsp.git
pip install git+ssh://git@github.com/MLGlobalHealth/sps.git
conda install conda-forge::arviz
conda install hydra-core
conda install jupyter
```

To test installations, type in terminal `python` and then `import numpyro`, `import jax`, `import sps`, `import dge`. Make sure none of these commands ends in an error.
