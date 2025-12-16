# About this repository

This repository contains codes for running [IsoQuant](https://github.com/ablab/IsoQuant) on RNA BaseCode bam file, and organise resulting files for an easier navigation/ visualisation.

## Run IsoQuant on RNA BaseCode bam file

Run IsoQuant on RNA BaseCode bam file through [docker image](https://hub.docker.com/r/basicgenomics/basecode_isoquant)

```bash
docker run --mount type=bind,src=$(pwd)/results/,dst=/usr/local/BaseCode/results/ --mount type=bind,src=$(pwd)/config/,dst=/usr/local/BaseCode/config/ --mount type=bind,src=/mnt/nvme/BaseCode_resources/,dst=/usr/local/app/resources/ --mount type=bind,src=/mnt/nvme/trial/fastq/,dst=/usr/local/BaseCode/fastq/ basicgenomics/basecode_isoquant:latest -j 10 --verbose
```

### Visualising IsoQuant outputs

#### Installation

1. Install Anaconda/[Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install)

2. After install Conda/Miniconda, run: 

```bash
conda env create -f isoquant_viewer.yaml
conda activate isoquant_viewer
conda install --channel conda-forge --channel bioconda pybedtools
```

Optionally, if you prefer R, install R in addition to the installation above. 
```bash
conda install  r-base -c conda-forge -y

```

3. Install [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) [v2.31.0](https://github.com/arq5x/bedtools2/releases/tag/v2.31.0)


The example notebooks are located in ```/utilities/notebook``` [^note]. 

[^note]:
    The codes are organised under the assumption that the notebook are exceuted from ```/utilities/notebook```, while the source codes are located in ```/utilities```.


