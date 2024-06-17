# sc_script_gen

I'm tired of looking up exactly how to process the data coming out of single cell assays.  I'm tired of copy/pasting
the same slurm job file 50 times and making minor alterations, hoping that this is the one I wrote for cellranger 8.0.0
and not 7.1.0.  I'm tired of juggling where this sample sheet or library file is supposed to go.

I already put all the information that is necessary in the demultiplexing sample sheet, so why should I need to do
anything else.

Maybe this would be better as part of a Nextflow pipeline, but I don't want to deal with the nightmare of debugging
Groovy code. Maybe Snakemake? Similar problem.  So, this thing.

# Installation

You're gonna need at least Python 3.10 and this is going to install [polars](https://pola.rs/),
[typer](https://typer.tiangolo.com/), and [rich](https://rich.readthedocs.io/en/stable/introduction.html).

```
pip install git+https://github.com/milescsmith/sc_script_gen
```

Recommend using [pipx](https://pipx.pypa.io/stable/installation/); if it is installed, just replace `pip` with `pipx`

# Usage

Right now, creating scripts for processing data from 10x Genomics Single Cell Immune Profiling (i.e. 5' scRNA-seq),
Single Cell ATAC-seq, and [ASAP-seq](https://cite-seq.com/asapseq/) are covered (though, note that the ASAP-seq
scripts use salmon alevin and not kallisto bustools because I don't like kallisto :shrug:).

Also note that there are currently some default values that are set to be convienent for me; I plan on moving those to
a `defaults` file in a later version.

There are two commands: `asapseq` and `fiveprime`:


## `asapseq`

Use information from a bcl-convert samplesheet to create scripts to process asapseq data using cellranger and asap_o_matic/salmon alevin

**Usage**:

```console
$ asapseq [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `create-asap-o-matic-script`: Batch create slurm scripts to concatenate...
* `create-atac-count-script`: Batch create slurm scripts to generate...
* `create-salmon-count-script`: Batch create slurm scripts to count...
* `version`: Prints the version of the package.

### `asapseq create-asap-o-matic-script`

Batch create slurm scripts to concatenate and rearrange ASAP-seq sequencing
data using asap-o-matic.

**Usage**:

```console
$ asapseq create-asap-o-matic-script [OPTIONS] SAMPLESHEET SCRIPTS_OUT_FOLDER FASTQ_PATH
```

**Arguments**:

* `SAMPLESHEET`: [required]
* `SCRIPTS_OUT_FOLDER`: [required]
* `FASTQ_PATH`: [required]

**Options**:

* `--demuxer [mkfastq|bcl-convert]`: [default: bcl-convert]
* `--conjugation [TotalSeqA|TotalSeqB]`: [default: TotalSeqB]
* `--num-cores INTEGER`: [default: 8]
* `--r2-reverse-complement / --no-r2-reverse-complement`: [default: no-r2-reverse-complement]
* `--mem INTEGER`: [default: 32]
* `--cpus INTEGER`: [default: 1]
* `--version`: Print version number.
* `--help`: Show this message and exit.

### `asapseq create-atac-count-script`

Batch create slurm scripts to generate count matrices for 10x Genomics scATAC-seq data

**Usage**:

```console
$ asapseq create-atac-count-script [OPTIONS] SAMPLESHEET SCRIPTS_OUT_FOLDER FASTQ_PATH
```

**Arguments**:

* `SAMPLESHEET`: [required]
* `SCRIPTS_OUT_FOLDER`: [required]
* `FASTQ_PATH`: [required]

**Options**:

* `--mem INTEGER`: [default: 32]
* `--cpus INTEGER`: [default: 1]
* `--ref PATH`: [default: /Volumes/shared-refs/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/]
* `--interval INTEGER`: [default: 5]
* `--max_jobs INTEGER`: [default: 8]
* `--mem_per_core INTEGER`: [default: 8]
* `-t, --template TEXT`: Path to the slurm cluster template  [default: /Volumes/guth_aci_informatics/software/slurm.template]
* `-l, --load_module`: Does the cellranger-atac module need to be loaded?  [default: True]
* `--version`: Print version number.
* `--help`: Show this message and exit.

### `asapseq create-salmon-count-script`

Batch create slurm scripts to count ASAP-seq protein data with salmon alevin

**Usage**:

```console
$ asapseq create-salmon-count-script [OPTIONS] SAMPLESHEET SCRIPTS_OUT_FOLDER RESULTS_PATH FASTQ_PATH
```

**Arguments**:

* `SAMPLESHEET`: Path to the bcl-convert samplesheet  [required]
* `SCRIPTS_OUT_FOLDER`: Path to where the generated scripts should be written  [required]
* `RESULTS_PATH`: Path to where the count results should be saved  [required]
* `FASTQ_PATH`: Path to the rearranged FASTQs  [required]

**Options**:

* `-i, --index PATH`: Path to the salmon index of the CITE-seq barcodes
* `-m, --mem INTEGER`: Amount of memory to use for each slurm job  [default: 32]
* `-c, --cpus INTEGER`: Amount of cpus to use for each slurm job  [default: 8]
* `--version`: Print version number.
* `--help`: Show this message and exit.

## `fiveprime`

**Usage**:

```console
$ fiveprime [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--help`: Show this message and exit.

**Commands**:

* `create_scripts`
* `version`: Prints the version of the package.

### `fiveprime create_scripts`

**Usage**:

```console
$ fiveprime create_scripts [OPTIONS] SAMPLESHEET SCRIPTS_OUT_FOLDER FASTQ_PATH
```

**Arguments**:

* `SAMPLESHEET`: Path to the bcl-convert samplesheet  [required]
* `SCRIPTS_OUT_FOLDER`: Path to where the generated scripts should be written  [required]
* `FASTQ_PATH`: Path to where the FASTQs are stored  [required]

**Options**:

* `--gex_index PATH`: Path to the salmon index of the CITE-seq barcodes  [default: /Volumes/shared-refs/cellranger/refdata-gex-GRCh38-2024-A]
* `--vdj_index PATH`: Path to the salmon index of the CITE-seq barcodes  [default: /Volumes/shared-refs/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0]
* `--feat_ref PATH`: Path to the salmon index of the CITE-seq barcodes  [default: /Volumes/guth_aci_informatics/references/miscellaneous/TotalSeq_C_Human_Universal_Cocktail_399905_Antibody_reference_UMI_counting.csv]
* `--chem [threeprime|SC3Pv1|SC3Pv2|SC3Pv3|SC3Pv4|SC3Pv3HT|SC-FB|fiveprime|SC5P-PE|SC5P-PE-v3|SC5P-R2|SC5P-R2-v3|SC5PHT|SFRP|MFRP|MFRP-R1|MFRP-RNA|MFRP-Ab|MFRP-Ab-R2pos50|MFRP-RNA-R1|MFRP-Ab-R1|ARC-v1]`: 10x kit chemistry. Currently, only 'SC5PHT' is supported  [default: SC5PHT]
* `--interval INTEGER`: [default: 2000]
* `--max_jobs INTEGER`: [default: 8]
* `--mem_per_core INTEGER`: [default: 8]
* `-m, --mem INTEGER`: Amount of memory to use for each slurm job  [default: 64]
* `-c, --cpus INTEGER`: Amount of cpus to use for each slurm job  [default: 8]
* `-t, --template TEXT`: Path to the slurm cluster template  [default: /Volumes/guth_aci_informatics/software/slurm.template]
* `-l, --load_module`: Does the cellranger-atac module need to be loaded?  [default: True]
* `--version`: Print version number.
* `--help`: Show this message and exit.