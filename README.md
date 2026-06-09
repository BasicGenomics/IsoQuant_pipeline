# BaseCode IsoQuant

Runs [IsoQuant](https://github.com/BasicGenomics/IsoQuant) (the Basic Genomics `basecode` fork, `v3.13.0.bg`) on the stitched‑molecule BAM produced by the [BaseCode Processing Pipeline](https://github.com/BasicGenomics/BaseCode). Utilities for visualising the IsoQuant output are included.

The fork adds a `--basecode` mode that makes IsoQuant compatible with BaseCode outputs. Without `--basecode` the tool is identical to upstream IsoQuant v3.13.0.

## Input

The pipeline expects the stitched‑molecule BAM produced upstream at:
```
results/<name>.stitched.molecules.sorted.bam
```
where `<name>` is the `name` set in `config/config.yaml`.

## Configuration

Edit `config/config.yaml`. Required:
- `name` — run name (used as the output prefix and `results/isoquant/<name>/` subfolder)
- `reference` — reference folder under `<resource_dir>/genome_references/<reference>/`

## Running

### Locally

Run from the repository root (relative paths resolve from there):

```bash
./BaseCodeIsoQuant -c config/config.yaml -t 30
```

Options:

| Flag | Description | Default |
| --- | --- | --- |
| `-c`, `--config` | Path to the config .yaml file | `config/config.yaml` |
| `-t`, `--threads` | Threads to allocate | `30` |
| `--verbose` | Stream full Snakemake / IsoQuant output | off |
| `--dry-run` | Print the execution plan without running | off |

Any unrecognised flags are forwarded straight to Snakemake.

Outputs are written to `results/isoquant/<name>/`; the MuData object is `results/isoquant/mudata/<name>_counts.h5mu`. Per‑rule logs are in `results/isoquant/logs/`.

### Docker

Available as a [Docker image](https://hub.docker.com/r/basicgenomics/basecode_isoquant). `cd` into the BaseCode run folder (where `results/` lives), then replace `/path/to/BaseCode_resources/` with your reference folder:

```bash
docker run \
  --mount type=bind,src=$(pwd)/results/,dst=/usr/local/BaseCodeIsoQuant/results/ \
  --mount type=bind,src=$(pwd)/config/,dst=/usr/local/BaseCodeIsoQuant/config/ \
  --mount type=bind,src=/path/to/BaseCode_resources/,dst=/usr/local/app/resources/ \
  --mount type=bind,src=$(pwd)/fastq/,dst=/usr/local/BaseCodeIsoQuant/fastq/ \
  basicgenomics/basecode_isoquant:latest -j 10 --verbose
```

## Visualising IsoQuant outputs

1. Install [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install) if you haven't already.

2. Create the environment:

   ```bash
   conda env create -f utilities/isoquant_viewer.yaml
   conda activate isoquant_viewer
   ```

   Optionally, for the R-based plots:

   ```bash
   conda install r-base -c conda-forge -y
   ```

Example notebooks are in `utilities/notebook/` and are meant to be run from that directory (the source code lives in `utilities/`).
