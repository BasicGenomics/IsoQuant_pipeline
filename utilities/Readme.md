## Visualising IsoQuant outputs

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

3. For viewing the genomic tracks, install [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) [v2.31.0](https://github.com/arq5x/bedtools2/releases/tag/v2.31.0) 


###

The example notebooks are located in ```/utilities/notebook``` [^note]. 

[^note]:
    The codes for visualising the results are organised under the assumption that the notebook are exceuted from ```/utilities/notebook```, while the source codes are located in ```/utilities```.


