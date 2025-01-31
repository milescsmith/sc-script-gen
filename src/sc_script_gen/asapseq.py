from enum import Enum
from pathlib import Path
from typing import Annotated, Optional

import typer
from polars.exceptions import ColumnNotFoundError
from rich import print as rprint
from rich.progress import Progress

from sc_script_gen import __version__
from sc_script_gen.utils import JOB_MANAGER_TEMPLATE_PATH, create_script_header, read_samplesheet

# obviously, this only works for me
# TODO: move this to a "defaults.toml" file?
CITESEQ_INDEX_PATH = "/Volumes/guth_aci_informatics/references/miscellaneous/salmon_totalseq_b_asapseq"
GENOMIC_REFERENCE_PATH = "/Volumes/shared-refs/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/"

# TODO: More feedback


class Conjugation(str, Enum):
    TotalSeqA = "TotalSeqA"
    TotalSeqB = "TotalSeqB"


class Demuxer(str, Enum):
    mkfastq = "mkfastq"
    bcl_convert = "bcl-convert"


asapseq = typer.Typer(
    name="ASAPseq script generator",
    help=(
        "Use information from a bcl-convert samplesheet to create scripts to process asapseq data using cellranger "
        "and asap_o_matic/salmon alevin"
    ),
    add_completion=False,
    no_args_is_help=True,
    add_help_option=True,
)


def version_callback(version: Annotated[bool, typer.Option("--version")] = True) -> None:  # FBT001
    """Prints the version of the package."""
    if version:
        rprint(f"[yellow]asapseq-script-gen[/] version: [bold blue]{__version__}[/]")
        raise typer.Exit()


def create_atac_count_script_body(
    sample: str,
    fastq_path: Path,
    reference_path: Path | None = None,
    job_interval: int = 2000,
    max_num_jobs: int = 8,
    mem_per_core: int = 8,
    job_template_path: str | None = None,
):
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
        f"asap-o-matic \\\n\t"
        f"--fastqs {fastq_path} \\\n\t"
        f"--sample {sample_id or f'{sample}_prot'} \\\n\t"
        f"--id {sample_id or f'{sample}_prot'} \\\n\t"
        f"--fastq_source {demuxer} \\\n\t"
        f"--conjugation {conjugation} \\\n\t"
        f"--outdir {fastq_path.joinpath(sample)} \\\n\t"
        f"--cores {num_cores}"
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
        index = CITESEQ_INDEX_PATH
    return (
        f"salmon alevin \\\n\t"
        f"--libType ISR \\\n\t"
        f"--index {index} \\\n\t"
        f"--mates1 {fastq_path.joinpath(sample, f'{sample_id}_R1.fastq.gz')} \\\n\t"
        f"--mates2 {fastq_path.joinpath(sample, f'{sample_id}_R2.fastq.gz')} \\\n\t"
        f"--output {results} \\\n\t"
        f"--threads {num_cores} \\\n\t"
        f"--citeseq \\\n\t"
        f"--featureStart 0 \\\n\t"
        f"--featureLength 15"
    )


# @app.callback(invoke_without_command=True)
@asapseq.command(no_args_is_help=True, context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
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
    reference_path: Annotated[Optional[Path], typer.Option("--ref")] = GENOMIC_REFERENCE_PATH,
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
    ] = JOB_MANAGER_TEMPLATE_PATH,
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
            is_eager=True,
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
        rprint(f"Found information for {len(samples)} samples")

    if reference_path is None:
        reference_path = GENOMIC_REFERENCE_PATH

    module_line = "module load cellranger-atac/2.1.0\n" if load_cellranger_module else "\n"

    with Progress() as progress_bar:
        task = progress_bar.add_task("Writing script files...", total=len(samples))
        for i in samples:
            progress_bar.console.print(f"Creating script for {i}")
            sample_fastq_path = fastq_path.joinpath(i)

            part1 = create_script_header(sample=i, jobtype="count", mem=mem, cpus=cpus)
            part2 = create_atac_count_script_body(
                sample=i,
                fastq_path=sample_fastq_path,
                reference_path=reference_path,
                job_interval=job_interval,
                max_num_jobs=max_num_jobs,
                mem_per_core=mem_per_core,
                job_template_path=job_template_path,
            )
            script_text = part1 + module_line + part2

            outfile = scripts_out_folder.joinpath(f"{i}_atac_count.job")
            progress_bar.console.print(f"Wrote script file to {outfile}")
            if not scripts_out_folder.resolve().exists():
                scripts_out_folder.resolve().mkdir()
            with outfile.open("w") as sf:
                sf.writelines(script_text)
            progress_bar.advance(task)


@asapseq.command(no_args_is_help=True, context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
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
        typer.Option("--version", callback=version_callback, help="Print version number.", is_eager=True),
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
    else:
        nsamples = ss_df[["Sample_ID", "Sample_Project"]].n_unique()
        rprint(f"Found information for {ss_df[['Sample_ID', 'Sample_Project']].n_unique()} samples")

    with Progress() as progress_bar:
        task = progress_bar.add_task("Writing script files...", total=nsamples)
        for i in ss_df[["Sample_ID", "Sample_Project"]].unique().iter_rows(named=True):
            progress_bar.console.print(f"Creating script for {i['Sample_ID']}")
            sample_fastq_path = fastq_path.joinpath(i["Sample_ID"])
            part1 = create_script_header(sample=i["Sample_ID"], jobtype="count", mem=mem, cpus=cpus)
            part2 = create_kite_script_body(
                sample=i["Sample_Project"],
                fastq_path=sample_fastq_path,
                sample_id=i["Sample_ID"],
                demuxer=demuxer.value,
                conjugation=conjugation.value,
                num_cores=num_cores,
                r2_reverse_complement=r2_reverse_complement,
            )
            script_text = part1 + part2
            outfile = scripts_out_folder.joinpath(f"{i['Sample_Project']}_asap-o-matic.job")
            progress_bar.console.print(f"Wrote script file to {outfile}")
            with outfile.open("w") as sf:
                sf.writelines(script_text)
            progress_bar.advance(task)


@asapseq.command(no_args_is_help=True)
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
        typer.Option("--version", callback=version_callback, help="Print version number.", is_eager=True),
    ] = False,
):
    """Batch create slurm scripts to count ASAP-seq protein data with salmon alevin"""
    # TODO: we should check to see if results_path exists and, if it doesn't, create it

    ss_df = read_samplesheet(samplesheet)
    ss_df = ss_df.filter(ss_df["Sample_ID"].str.ends_with("prot"))

    if "Sample_Project" not in ss_df.columns:
        msg = "'Sample_Project' was not found in the samplesheet columns. Is your samplesheet in the correct format?"
        raise ColumnNotFoundError(msg)
    else:
        nsamples = ss_df[["Sample_ID", "Sample_Project"]].n_unique()
        rprint(f"Found information for {nsamples} samples")

    # seems like this should be extracted into a new method
    # and those progress_bar parts should be a decorator?
    with Progress() as progress_bar:
        task = progress_bar.add_task("Writing script files...", total=nsamples)
        for i in ss_df[["Sample_ID", "Sample_Project"]].unique().iter_rows(named=True):
            progress_bar.console.print(f"Creating script for {i['Sample_ID']}")
            part1 = create_script_header(sample=i["Sample_ID"], jobtype="alevin", mem=mem, cpus=cpus)
            part2 = create_salmon_script_body(
                sample=i["Sample_Project"],
                fastq_path=fastq_path,
                sample_id=i["Sample_ID"],
                results=results_path.joinpath(i["Sample_Project"]),
                index=index_path,
                num_cores=cpus,
            )

            script_text = part1 + part2
            outfile = scripts_out_folder.joinpath(f"{i['Sample_Project']}_salmon_count.job")
            progress_bar.console.print(f"Wrote script file to {outfile}")
            with outfile.open("w") as sf:
                sf.writelines(script_text)
            progress_bar.advance(task)
