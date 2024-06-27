# VAE_PopGen

## Conda environment

```
conda create --name popgen python=3.12
conda activate popgen
pip install "jax[cpu]"
conda install conda-forge::numpyro
conda install conda-forge::arviz
conda install hydra-core
pip install git+ssh://git@github.com/MLGlobalHealth/dge.git
pip install git+ssh://git@github.com/MLGlobalHealth/sps.git
```

If working on a machine with CPU, replace the lines installing `jax` and `numpyro` with 

```
pip install --upgrade "jax[cuda12]"
pip install 'numpyro[cuda]' -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
```

To test installations, type in terminal `python` and then `import numpyro`, `import jax`, `import sps`, `import dge`. Make sure none of these commands ends in an error.