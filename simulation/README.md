## 1. Environment Configuration

### 1.1 Install dependencies

```shell
# Create enviroment
mamba create -n simulation
mamba activate simulation

# VISOR
mamba install VISOR

# bcftools
mamba install bcftools

# samtools
mamba install samtools

# scipy
mamba install scipy

# pysam
mamba install pysam

# cython
pip install cython

# TEMP3 is also needed
```

### 1.2 Compile cython scripts

```shell
# Current working directory should be your_path_to/simulation
python setup.py build_ext -i && rm -r build && rm *c
```

## 2. Simulation with blades

To improve efficiency, the simulation steps can be performed vi blades, following `simulate_training_data` or `simulate_test_data`.