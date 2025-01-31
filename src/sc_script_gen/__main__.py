import typer

from sc_script_gen.asapseq import asapseq
from sc_script_gen.scrnaseq import scrnaseq

script_gen = typer.Typer(
    name="10X processing script generator",
    help=(
        "Create the scripts necessary to process raw sequencing data coming from 10X Genomics assays "
        "into count matrices"
    ),
    add_completion=False,
    no_args_is_help=True,
    add_help_option=True,
)
script_gen.add_typer(
    typer_instance=asapseq,
    name="asapseq",
    no_args_is_help=True,
    help=(
        "Use information from a bcl-convert samplesheet to create scripts to process asapseq data using cellranger "
        "and asap_o_matic/salmon alevin"
    ),
)
script_gen.add_typer(
    typer_instance=scrnaseq,
    name="scrnaseq",
    no_args_is_help=True,
    help="Use information from a bcl-convert samplesheet to create scripts to process data from the 10x Genomics Single Cell Immune Profiling data using cellranger",
)
