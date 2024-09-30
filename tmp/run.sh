module load mamba
mamba env create --file devtools/conda-envs/test-env.yaml
module load gsl-2.7.1-gcc-12.1.0
pip install -e .

python tmp/run.py