# Standard LDSC environment

The primary comparison runs **original bulik/ldsc 1.0.1**, commit
`2fdeeb3b44379408794154993dbd6101b8946b7e`, directly under **Python 2.7.18**.
Neither `ldsc.py`, `munge_sumstats.py`, nor their statistical modules are patched
or translated. Keep this environment separate from GCMSEM's Python 3 environment.

The executed numerical dependencies are NumPy 1.16.0, SciPy 0.18.0, pandas
0.20.3 and bitarray 0.8.3. They satisfy the original repository's dependency
ranges. [standard-ldsc-requirements.txt](standard-ldsc-requirements.txt) pins the
complete regression runtime. `pybedtools`, used by upstream `make_annot.py`, is
unnecessary here because this tutorial uses published LD scores.

## Upstream Conda route

From the GCMSEM repository root, fetch the untouched source and use the
[upstream installation procedure](https://github.com/bulik/ldsc#readme):

```bash
python tutorials/setup_standard_ldsc.py --out external/standard-ldsc --fetch-only
conda env create --file external/standard-ldsc/environment.yml
conda activate ldsc
python -m pip install pip==20.3.4 setuptools==44.1.1 wheel==0.34.2
python -m pip install -r tutorials/standard-ldsc-requirements.txt
gcmsem_ldsc_python="$(python -c 'import sys; print(sys.executable)')"
conda deactivate
python tutorials/setup_standard_ldsc.py --out external/standard-ldsc \
  --python2 "$gcmsem_ldsc_python"
```

The upstream environment is historical: it needs compatible legacy packages
and a compatible platform. These exact versions were executed on **Linux
x86-64**; availability on other platforms is not established. If the historical
Conda solve is unavailable, the following source-build route reproduces the
runtime used for this release.

## Tested Python source-build route

On Ubuntu/Linux x86-64, install the compiler and development headers first:

```bash
sudo apt-get install build-essential libbz2-dev zlib1g-dev libreadline-dev libffi-dev libssl-dev
```

The release was executed with unmodified CPython 2.7.18 from python.org, GCC
13.3.0, official PyPI numerical wheels, and a locally compiled bitarray
extension. The Python source archive has SHA-256
`da3080e3b488f648a3d7a4560ddee895284c3380b11d6de75edb986526b9a814`.

Run from the GCMSEM repository root:

```bash
gcmsem_python_build="$PWD/external/python-build"
gcmsem_python_prefix="$PWD/external/python2"
mkdir -p "$gcmsem_python_build"
curl -fL https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz \
  -o "$gcmsem_python_build/Python-2.7.18.tgz.part"
printf '%s  %s\n' da3080e3b488f648a3d7a4560ddee895284c3380b11d6de75edb986526b9a814 \
  "$gcmsem_python_build/Python-2.7.18.tgz.part" | sha256sum --check - || exit 1
mv "$gcmsem_python_build/Python-2.7.18.tgz.part" "$gcmsem_python_build/Python-2.7.18.tgz"
tar --no-same-owner -xzf "$gcmsem_python_build/Python-2.7.18.tgz" -C "$gcmsem_python_build"
(
  cd "$gcmsem_python_build/Python-2.7.18"
  ./configure --prefix="$gcmsem_python_prefix" --enable-unicode=ucs4 --with-ensurepip=install
  make -j4
  make install
)
gcmsem_ldsc_python="$gcmsem_python_prefix/bin/python2.7"
"$gcmsem_ldsc_python" -m pip install pip==20.3.4 setuptools==44.1.1 wheel==0.34.2
"$gcmsem_ldsc_python" -m pip install -r tutorials/standard-ldsc-requirements.txt
python tutorials/setup_standard_ldsc.py --out external/standard-ldsc \
  --python2 "$gcmsem_ldsc_python"
```

Install `libbz2-dev` before building Python: stock `munge_sumstats.py` imports
`bz2`. In the execution environment the header was added after the first
build, so the **unmodified CPython bzip2 extension** was compiled separately;
this did not change LDSC code. Numerical wheels were downloaded separately and
installed offline there; the commands above install the same pinned packages.

The setup helper verifies the pinned commit, rejects tracked source edits,
checks both stock commands, and writes source hashes and dependency versions
to `external/standard-ldsc-setup.json`. The executed environment passed **125
upstream tests** covering regression, IRWLS, jackknife, parsing and summary
statistics. Pass `--python2 "$gcmsem_ldsc_python"` to the tutorial runner.

## Munging chunk size

The runner supplies stock `munge_sumstats.py --chunksize 500000` to bound input
parsing batches. This avoids a slow array-membership path in the historical
pandas/NumPy combination. It does not change the QC thresholds or regression
settings. The source, allele checks, published reference/weight files and
LDSC statistical defaults remain unchanged. The stock `ldsc.py --rg` command
receives no custom statistical options.

As an independent check, the complete GSCAN drinks/week input was munged with
both 250,000 and 500,000 rows per chunk. The decompressed outputs were
**byte-identical**: SHA-256
`760214a0820658e99ffa898c5fa6217ae64ba7ebec840d117091c0af2f61f6e0`.
Compressed bytes can differ because stock gzip headers contain timestamps.
