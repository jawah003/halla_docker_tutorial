# Installing HAllA natively on macOS

## The recipe

Everything goes into a dedicated conda environment. Do not install into `base`.

### 1. Environment with Python and R together

```bash
conda create -n halla_native python=3.9 pip -y
conda activate halla_native
conda install -c conda-forge r-base=4.3 r-eva r-xicor -y
```


### 2. HAllA, then rpy2 rebuilt from source

```bash
pip install halla==0.8.40
```

Then force `rpy2` to compile against the conda R rather than using a prebuilt wheel:

```bash
LDFLAGS="-L$CONDA_PREFIX/lib/R/lib" \
CPPFLAGS="-I$CONDA_PREFIX/lib/R/include" \
pip install rpy2 rpy2-rinterface --no-binary :all: --no-cache-dir --force-reinstall
```

### 3. Pin the R paths to the environment

Without this, the variables drift back to system defaults the next time you open a shell. These
scripts set them on activation and clear them on deactivation:

```bash
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d" "$CONDA_PREFIX/etc/conda/deactivate.d"

cat > "$CONDA_PREFIX/etc/conda/activate.d/r_paths.sh" <<'INNER'
export R_HOME="$CONDA_PREFIX/lib/R"
export DYLD_LIBRARY_PATH="$CONDA_PREFIX/lib/R/lib:$DYLD_LIBRARY_PATH"
INNER

cat > "$CONDA_PREFIX/etc/conda/deactivate.d/r_paths.sh" <<'INNER'
unset R_HOME
unset DYLD_LIBRARY_PATH
INNER
```

### 4. Pin Pillow

```bash
pip install "Pillow<10.0.0"
```

Pillow 10 removed a deprecated `mode` argument that the matplotlib in HAllA's dependency set still
passes. Without the pin, a run computes everything correctly and then crashes while writing the
hallagram PDF.

### 5. Verify

Reactivate so the path scripts run, then check both halves:

```bash
conda deactivate && conda activate halla_native

halla -h
python -c "import rpy2.robjects.packages as r; print('eva:', r.isinstalled('eva')); print('XICOR:', r.isinstalled('XICOR'))"
```

`halla -h` should print the usage block with no R linkage errors, and both packages should report
`True`. If `halla -h` works but the packages report `False`, `rpy2` has loaded a different R than
the one conda installed — check `R_HOME` and redo step 3.

## Running the tutorial without Docker

The notebook shells out to `docker run`. To use a native install instead, replace the command in
`build_halla_cmd` with the bare CLI call:

```python
def build_halla_cmd(name, cfg):
    cmd = ["halla", "-x", "data/X_metabolites.txt", "-y", "data/Y_taxa.txt", ...]
```

Drop the `docker run --rm -v ... -w /halla IMAGE` prefix and keep everything from `"halla"` onward.
Paths are then relative to your working directory rather than to the container mount, so run the
notebook from the repository root. Skip Section 3 (the image build) entirely.
