
# BE3D (BEClust3D)

**BE3D (or BEClust3D)** is a Python package for interpreting structure-function relationships in base editor (BE) tiling mutagenesis data. The workflow includes 3 main modules: (1) quality assessment and statistical analysis of screen data by gene, (2) extrapolation of BE screen signals onto 3D structures, identification of significant residues, and clustering to identify hotspots from a structure-function perspective, and (3) aggregation of multiple screens for the idenficiation of signficiant residues and clusters

You can run the BE3D pipeline either:
- On **Google Colab** (no installation required), or
- **Locally** (faster execution)

---

## Overview

- **[Workflow Overview](#workflow-overview)**
- **[Input](#input)**
- **[Features](#features)**
    - **[Quality Assessment](#quality-assessment-hypothesis-test-visualization)**
    - **[3D-based Scoring and Visualization](#be-clust3d-visualization-of-lfc-and-lfc3d-hits)**
    - **[Meta-Aggregation](#be-metaclust3d)**
    - **[Genomics 2 Portal Output](#visualization-on-the-genomics-2-portal)**
- **[Installation](#installation)**
- **[Examples](#getting-started-examples)**
- **[Some Extra Notes](#notes)**
- **[Google Colab Notebooks](#sample-google-colab-notebooks)**
- **[License](#license)**

## Workflow Overview

The following figure provides an overview of the BE3D workflow:

![BE3D workflow](imgs/BE3D_workflow.png)

BE3D enables structure-function analysis of BE tiling mutagenesis data by mapping mutation readouts (log fold change, LFC) onto 3D protein structures. This can be extended to multiple screens or cross-species comparisons. The workflow consists of:

**A. BE-QA**: Assesses the quality of BE screens by testing if annotated knockout guides (e.g., nonsense or splice site) and annotated neutral guides (e.g., silent or no mutation) guides have significantly different LFC score distributions.

**B. BE-Clust3D**: Maps LFC values by amino acid residue onto 3D protein structures and computes a per-residue 3D-normalized LFC score (LFC3D score) based on spatial proximity (default: 6 Å). Then, agglomerative clustering is performed with a second spatial proximity parameter (default: 6 Å) to identify hotspots of potential functional importance. This clustering can be performed on either the original LFC values by amino acid, or the new LFC3D score. 

**C. BE-MetaClust3D**: Aggregates data from multiple screens to enhance signal strength and detect residues that might be missed due to the noise present in BE screens. Cross-species or cross-isoform screens can also be integrated together through an optional alignment step. 

## Input

BE3D requires the following inputs:

1. **BE Screen Scores (TSV)**: Must include Mutation Category, Amino Acid Edit, Gene Name, and Score. You must indicate column names as part of the input. 

    Example TSV:

    ```tsv
    Gene	sgRNA_score	mutation	predicted_edits
    MEN1	-0.18977	Missense	Gly2Arg;Met1Ile
    MEN1	-0.22247	Silent	Leu10Leu
    ```

    Example input config (Python):

    ```python
    mut_col   = "mutation"
    val_col   = "sgRNA_score"
    gene_col  = "Gene"
    edits_col = "predicted_edits"
    ```

2. **Uniprot ID**: Required to fetch protein sequence and structure from UniProt/AlphaFold.

    ```python
    input_uniprot = "O00255" # (MEN1)
    ```

4. **Optional FASTA and PDB**: Provide custom protein sequence and structure files for non-canonical sequences, non-canonical proteins, or alternative structures. If these fields are left empty, the pipeline fetches the AlphaFold of the canonical isoform structure for the given Uniprot ID. 

    ```python
    input_pdb   = 'men1_AF3.pdb'
    input_fasta = 'men1.fasta'
    ```

## Features

### Quality Assessment: Hypothesis Test Visualization

BE-QA performs Mann-Whitney and Kolmogorov-Smirnov tests on LFC distributions, comparing knockout and neutral mutations. Knockout mutations of a single gene in a single screen are compared against neutral mutations of that single gene (hypothesis 1) or neutral mutations of all genes in that screen (hypotehsis 2). Results are visualized with statistical annotations.

![QA](imgs/QA.png)

### BE-Clust3D: Visualization of LFC and LFC3D Hits

BE-Clust3D prioritizes residues by aggregating LFC values within a defined spatial range. This enhances signal detection by extrapolating functional data in 3D space to calculate an LFC3D score. Results are visualized and clustered. 

This step also includes the preprocessing of scores organized by sgRNA to scores organized by residues, running sequence alignment to combine screens on different genes, and the calculation of p-values to define statistical thresholds to define what is a hit. 

![LFC/LFC3D](imgs/LFC_and_LFC3D.png)

### BE-MetaClust3D

BE-MetaClust3D aggregates across multiple screens to identify consensus hotspots and enhance weaker signals across multiple screens. 

![Meta-Aggregation](imgs/Meta-aggregation.png)

### Visualization on the Genomics 2 Portal

Results are provided in G2P-compatible TSV file, which can be downloaded and interactively viewable in [Genomics 2 Portal](https://g2p.broadinstitute.org/mapping).

![G2P](imgs/G2P.png)

## Installation
### 1. Install BE3D

Install directly from GitHub:

``` bash
pip install git+https://github.com/broadinstitute/beclust3d-public.git
```

> This installation command is also included in the example Google Colab
> notebooks.

### 2. Create Conda Environment

We recommend creating a dedicated conda environment before running BE3D
locally.

#### Apple Silicon (ARM Mac OS)

``` bash
conda env create -f environment_arm.yml
```

#### Linux (x86_64)

``` bash
conda env create -f environment.yml
```

Activate the environment:

``` bash
conda activate beclust3d
```

## Getting Started
### Running BE3D Locally

The script `examples/beclust3d_local.py` runs BE3D using a YAML
configuration file that specifies:

-   Input screen data
-   Structural model
-   Parameters
-   Output directory

Example usage:

``` bash
conda activate beclust3d
cd examples/

# DNMT3A example (Lue et al.)
python beclust3d_local.py yaml/dnmt3a_local.yaml

# MEN1 example (Perner et al.)
python beclust3d_local.py yaml/men1_local.yaml
```

## Google Colab Notebooks
BE3D can also be run directly in Google Colab.

**Single Screen Notebook Example (DNMT3A)**: 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/broadinstitute/BE3D/blob/main/examples/BE3Dv6-SingleScreen-Notebook.ipynb)

**Multi Screen Notebook with Meta-Aggregation Example (MEN1)**: 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/broadinstitute/BE3D/blob/main/examples/BE3Dv6-MultipleScreens-Notebook.ipynb)

**Multi Screen Notebook with Meta-Aggregation and Conservation Example**: 
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/broadinstitute/BE3D/blob/main/examples/BE3Dv6-MultipleScreensConservation-Notebook.ipynb)

## Github Structure

## Notes

### Structure

The pipeline automatically queries the UNIPROT protein sequence and AlphaFold structure of the protein of interest. If users want to use a PDB or other custom structure, they would need to upload the ```structure.pdb``` file and provide the filepath to the structure. 

The pipeline also automatically uses DSSP to annotate a pdb file for secondary structures. However, this tool is known to sometimes fail on larger structures. For a custom PDB upload, it is recommended that the user uploads their own DSSP file, as DSSP may fail on these structures. The annotations for DSSP are not necessary for the pipeline until the final characterization step, and would not affect preprocessing, prioritizing hits, meta-aggregation, or clustering. Even for the final characterization step, DSSP annotations are an optional input. 

The DSSP Web Portal can be found at: https://pdb-redo.eu/dssp

### Conservation

For sequence alignment, the pipeline runs MUSCLE locally in order to align 2 sequences in order to compare between isoforms or across species. 

For running CLUSTAL, the associated formating packages do not work for arm machines (ie M1/M2/M3 MacBooks). However, the packages should download for Windows and Linux based machines. If the user is using an arm machine, it is recommended to set the mode to 'query' instead of 'run', which calls the MUSCLE API. 

If MUSCLE or CLUSTAL cannot be run locally, the pipeline queries the MUSCLE API, although this may also fail due to issues with the API. Running the MUSCLE API also skips the next step using CLUSTAL. 

Another option to skip MUSCLE and CLUSTAL is for users to run alignment on their own in a CLUSTAL format, and provide the ```sequence.align``` alignment file into the pipeline which is one of the optional inputs. 



## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE.txt) file for details.
