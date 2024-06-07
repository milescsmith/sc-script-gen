from enum import Enum
from pathlib import Path
from typing import Annotated, Optional

import polars as pl
import typer
from polars.exceptions import ColumnNotFoundError
from rich import print as rprint

from asapseq_script_gen import __version__

# obviously, this only works for me
# TODO: move this to a "defaults.toml" file?
JOB_MANAGER_TEMPLATE_PATH = "/Volumes/guth_aci_informatics/software/slurm.template"
CITESEQ_REFERENCE_PATH = "/Volumes/guth_aci_informatics/references/miscellaneous/salmon_totalseq_b_asapseq"
GENOMIC_REFERENCE_PATH = "/Volumes/shared-refs/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/"

# TODO: More feedback

class Conjugation(str, Enum):
    TotalSeqA = "TotalSeqA"
    TotalSeqB = "TotalSeqB"


class Demuxer(str, Enum):
    mkfastq = "mkfastq"
    bcl_convert = "bcl-convert"


app = typer.Typer(
    name="script gen",
    help="Use information from a bcl-convert samplesheet to create scripts to process asapseq data using cellranger and asap_o_matic/salmon alevin",
)

@app.command(name="version", context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def version_callback(value: Annotated[bool, typer.Option()] = True) -> None:  # FBT001
    """Prints the version of the package."""
    if value:
        rprint(f"[yellow]asapseq-script-gen[/] version: [bold blue]{__version__}[/]")
        raise typer.Exit()


def read_samplesheet(samplesheet_file: Path) -> pl.DataFrame:
    """This assumes we're using a bcl-convert style samplesheet"""
    with samplesheet_file.open("r") as s:
        n = 1
        while "[BCLConvert_Data]" not in s.readline():
            n += 1

    return pl.read_csv(samplesheet_file, skip_rows=n)


def create_script_header(
    sample: str,
    jobtype: str | None = None,
    mem: int = 32,
    cpus: int = 1,
):
    """Generate a standard slurm job file header"""
    return (
        f"#! /bin/bash -l\n\n"
        f"#SBATCH -J {sample}{'_' + jobtype if jobtype else ''}\n"
        f"#SBATCH -o {sample}{'_' + jobtype if jobtype else ''}.log\n"
        f"#SBATCH --mail-user=None\n"
        f"#SBATCH --mail-type=END,FAIL\n"
        f"#SBATCH --mem={mem}G\n"
        f"#SBATCH --partition=serial\n"
        f"#SBATCH --nodes=1\n"
        f"#SBATCH --cpus-per-task={cpus}\n"
    )


def create_atac_count_script_body(
    sample: str,
    fastq_path: Path,
    reference_path: Path | None = None,
    job_interval: int = 5,
    max_num_jobs: int = 8,
    mem_per_core: int = 8,
    job_template_path: str | None = None,
):
    if job_template_path is None:
        job_template_path = JOB_MANAGER_TEMPLATE_PATH
    if reference_path is None:
        reference_path = GENOMIC_REFERENCE_PATH

    return (
        f"cellranger-atac \\\n"
        f"\tcount \\\n"
        f"\t\t--id={sample} \\\n"
        f"\t\t--reference={reference_path} \\\n"
        f"\t\t--fastqs={fastq_path} \\\n"
        f"\t\t--jobmode={job_template_path} \\\n"
        f"\t\t--disable-ui \\\n"
        f"\t\t--jobinterval={job_interval} \\\n"
        f"\t\t--maxjobs={max_num_jobs} \\\n"
        f"\t\t--mempercore={mem_per_core} \\\n"
        f"\t\t--sample={sample}_scATAC"
    )


def create_kite_script_body(
    sample: str,
    fastq_path: Path,
    sample_id: str | None = None,
    demuxer: Demuxer = Demuxer.bcl_convert,
    conjugation: Conjugation = Conjugation.TotalSeqB,
    num_cores: int = 8,
    r2_reverse_complement: bool = True,
) -> str:
    asap_o_matic_script = (
        "asap-o-matic \\\n"
        f"\t--fastqs /s/guth-aci/ARA08/data/raw/fastqs/asapseq_set_2/{sample} \\\n"
        f"\t--sample {sample_id if sample_id else sample + '_prot'} \\\n"
        f"\t--id {sample_id if sample_id else sample + '_prot'} \\\n"
        f"\t--fastq_source {demuxer} \\\n"
        f"\t--conjugation {conjugation} \\\n"
        f"\t--outdir {fastq_path.joinpath(sample)} \\\n"
        f"\t--cores {num_cores}"
    )
    if not r2_reverse_complement:
        asap_o_matic_script += " \\\n\t--no-rc-R2"

    return asap_o_matic_script


def create_salmon_script_body(
    sample: str,
    fastq_path: Path,
    sample_id: str,
    results: Path,
    index: Path | None = None,
    num_cores: int = 8,
) -> str:
    if index is None:
        index = CITESEQ_REFERENCE_PATH
    return (
        "salmon alevin \\\n"
        "\t--libtype ISR \\\n"
        f"\t--index {index} \\\n"
        f"\t--mates1 {fastq_path.joinpath(sample, sample_id+'_R1.fastq.gz')} \\\n"  # /s/guth-aci/ana_multiome_asapseq/data/raw/fastqs/set1/scATAC_1/scATAC_1_prot_R1.fastq.gz \\\n"
        f"\t--mates2 {fastq_path.joinpath(sample, sample_id+'_R2.fastq.gz')} \\\n"  # /s/guth-aci/ana_multiome_asapseq/data/raw/fastqs/set1/scATAC_1/scATAC_1_prot_R2.fastq.gz \\\n"
        f"\t--output {results} \\\n"
        f"\t--threads {num_cores} \\\n"
        "\t--citeseq \\\n"  # these next three lines are invariant for CITE-seq. Only need to change this if we aren't doing ASAPseq in which case, why are you using this?
        "\t--featureStart 0 \\\n"
        "\t--featureLength 15"
    )


# @app.callback(invoke_without_command=True)
@app.command(no_args_is_help=True)
def create_atac_count_script(
    # Info used in the loop
    samplesheet: Annotated[
        Path,
        typer.Argument(),
    ],
    scripts_out_folder: Annotated[
        Path,
        typer.Argument(),
    ],
    fastq_path: Annotated[Path, typer.Argument()],
    # universal settings for the cellranger martian jobmanager
    mem: Annotated[
        int,
        typer.Option("--mem"),
    ] = 32,
    cpus: Annotated[
        int,
        typer.Option("--cpus"),
    ] = 1,
    reference_path: Annotated[Optional[Path], typer.Option("--ref")] = None,
    job_interval: Annotated[
        int,
        typer.Option("--interval"),
    ] = 5,
    max_num_jobs: Annotated[
        int,
        typer.Option("--max_jobs"),
    ] = 8,
    mem_per_core: Annotated[
        int,
        typer.Option("--mem_per_core"),
    ] = 8,
    job_template_path: Annotated[
        Optional[str],
        typer.Option("--template", "-t", help="Path to the slurm cluster template"),
    ] = None,
    load_cellranger_module: Annotated[
        bool,
        typer.Option(
            "--load_module",
            "-l",
            help="Does the cellranger-atac module need to be loaded?",
        ),
    ] = True,
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            callback=version_callback,
            help="Print version number.",
        ),
    ] = False,
) -> None:
    """Batch create slurm scripts to generate count matrices for 10x Genomics scATAC-seq data"""
    ss_df = read_samplesheet(samplesheet)

    if "Sample_Project" not in ss_df.columns:
        msg = "'Sample_Project' was not found in the samplesheet columns. Is your samplesheet in the correct format?"
        raise ColumnNotFoundError(msg)
    else:
        samples = ss_df["Sample_Project"].unique()

    if reference_path is None:
        reference_path = GENOMIC_REFERENCE_PATH

    module_line = "module load cellranger-atac/2.1.0\n" if load_cellranger_module else "\n"

    for i in samples:
        sample_fastq_path = fastq_path.joinpath(i)

        script_text = (
            f"{
                create_script_header(
                    sample=i,
                    jobtype='count',
                    mem=mem,
                    cpus=cpus
                )
            }"
            f"{module_line}"
            f"{
                create_atac_count_script_body(
                    sample=i,
                    fastq_path=sample_fastq_path,
                    reference_path=reference_path,
                    job_interval=job_interval,
                    max_num_jobs=max_num_jobs,
                    mem_per_core=mem_per_core,
                    job_template_path=job_template_path
                )
            }"
        )
        outfile = scripts_out_folder.joinpath(f"{i}_atac_count.job")
        with outfile.open("w") as sf:
            sf.writelines(script_text)


@app.command(no_args_is_help=True)
def create_asap_o_matic_script(
    samplesheet: Annotated[
        Path,
        typer.Argument(),
    ],
    scripts_out_folder: Annotated[
        Path,
        typer.Argument(),
    ],
    fastq_path: Annotated[Path, typer.Argument()],
    demuxer: Annotated[Demuxer, typer.Option()] = Demuxer.bcl_convert,
    conjugation: Annotated[Conjugation, typer.Option()] = Conjugation.TotalSeqB,
    num_cores: Annotated[int, typer.Option()] = 8,
    r2_reverse_complement: Annotated[bool, typer.Option()] = False,
    mem: Annotated[
        int,
        typer.Option("--mem"),
    ] = 32,
    cpus: Annotated[
        int,
        typer.Option("--cpus"),
    ] = 1,
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            callback=version_callback,
            help="Print version number.",
        ),
    ] = False,
):
    """Batch create slurm scripts to concatenate and rearrange ASAP-seq sequencing
    data using asap-o-matic.
    """
    ss_df = read_samplesheet(samplesheet)
    ss_df = ss_df.filter(ss_df["Sample_ID"].str.ends_with("prot"))

    if "Sample_Project" not in ss_df.columns:
        msg = "'Sample_Project' was not found in the samplesheet columns. Is your samplesheet in the correct format?"
        raise ColumnNotFoundError(msg)

    for i in ss_df[["Sample_ID", "Sample_Project"]].unique().iter_rows(named=True):
        sample_fastq_path = fastq_path.joinpath(i["Sample_ID"])
        script_text = (
            f"{
                create_script_header(
                    sample=i['Sample_ID'],
                    jobtype="count",
                    mem=mem,
                    cpus=cpus
                )
            }"
            f"{
                create_kite_script_body(
                    sample=i['Sample_Project'],
                    fastq_path=sample_fastq_path,
                    sample_id=i['Sample_ID'],
                    demuxer=demuxer.value,
                    conjugation=conjugation.value,
                    num_cores=num_cores,
                    r2_reverse_complement=r2_reverse_complement
                )
            }"
        )
        outfile = scripts_out_folder.joinpath(f"{i['Sample_Project']}_atac_count.job")
        with outfile.open("w") as sf:
            sf.writelines(script_text)


@app.command(no_args_is_help=True)
def create_salmon_count_script(
    samplesheet: Annotated[
        Path,
        typer.Argument(help="Path to the bcl-convert samplesheet"),
    ],
    scripts_out_folder: Annotated[
        Path,
        typer.Argument(help="Path to where the generated scripts should be written"),
    ],
    results_path: Annotated[Path, typer.Argument(help="Path to where the count results should be saved")],
    fastq_path: Annotated[Path, typer.Argument(help="Path to the rearranged FASTQs")],
    index_path: Annotated[
        Optional[Path],
        typer.Option("--index", "-i", help="Path to the salmon index of the CITE-seq barcodes"),
    ] = None,
    mem: Annotated[
        int,
        typer.Option("--mem", "-m", help="Amount of memory to use for each slurm job"),
    ] = 32,
    cpus: Annotated[
        int,
        typer.Option("--cpus", "-c", help="Amount of cpus to use for each slurm job"),
    ] = 8,
    version: Annotated[
        bool,
        typer.Option(
            "--version",
            callback=version_callback,
            help="Print version number.",
        ),
    ] = False,
):
    """Batch create slurm scripts to count ASAP-seq protein data with salmon alevin"""
    # TODO: we should check to see if results_path exists and, if it doesn't, create it

    ss_df = read_samplesheet(samplesheet)
    ss_df = ss_df.filter(ss_df["Sample_ID"].str.ends_with("prot"))

    if "Sample_Project" not in ss_df.columns:
        msg = "'Sample_Project' was not found in the samplesheet columns. Is your samplesheet in the correct format?"
        raise ColumnNotFoundError(msg)

    for i in ss_df[["Sample_ID", "Sample_Project"]].unique().iter_rows(named=True):
        script_text = (
            f"{
                create_script_header(
                    sample=i['Sample_ID'],
                    jobtype='alevin',
                    mem=mem,
                    cpus=cpus
                )
            }"
            f"{
                create_salmon_script_body(
                    sample=i['Sample_Project'],
                    fastq_path=fastq_path,
                    sample_id=i['Sample_ID'],
                    results=results_path,
                    index=index_path,
                    num_cores=cpus,
                )
            }"
        )
        outfile = scripts_out_folder.joinpath(f"{i['Sample_Project']}_salmon_count.job")
        with outfile.open("w") as sf:
            sf.writelines(script_text)
