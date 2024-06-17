from enum import Enum
from pathlib import Path
from typing import Annotated, Literal, Optional

import typer
from polars.exceptions import ColumnNotFoundError
from rich import print as rprint

from sc_script_gen import __version__
from sc_script_gen.utils import JOB_MANAGER_TEMPLATE_PATH, create_script_header, read_samplesheet

GEX_REFERENCE_PATH = Path("/Volumes/shared-refs/cellranger/refdata-gex-GRCh38-2024-A")
VDJ_REFERENCE_PATH = Path("/Volumes/shared-refs/cellranger/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0")
FEAT_REFERENCE_PATH = Path(
    "/Volumes/guth_aci_informatics/references/miscellaneous/TotalSeq_C_Human_Universal_Cocktail_399905_Antibody_reference_UMI_counting.csv"
)

# NOTE: right now, this is only going to cover the options that I actually use
# more options will probably not be covered until I 1) need them or 2) am sufficiently bored


class Chemistry(str, Enum):
    threeprime = "threeprime"
    SC3Pv1 = "SC3Pv1"  # Single Cell 3' v1, v2, v3, or v4
    SC3Pv2 = "SC3Pv2"
    SC3Pv3 = "SC3Pv3"
    SC3Pv4 = "SC3Pv4"
    SC3Pv3HT = "SC3Pv3HT"  # Single Cell 3' v3.1 HT
    SC_FB = "SC-FB"  # Single Cell Antibody-only 3' v2 or 5'
    fiveprime = "fiveprime"  # Single Cell 5'
    SC5P_PE = "SC5P-PE"  # Paired-end Single Cell 5'
    SC5P_PE_v3 = "SC5P-PE-v3"  # Paired-end Single Cell 5' v3
    SC5P_R2 = "SC5P-R2"  # R2-only Single Cell 5'
    SC5P_R2_v3 = "SC5P-R2-v3"  # R2-only Single Cell 5' v3
    SC5PHT = "SC5PHT"  # Single Cell 5' v2 HT
    SFRP = "SFRP"  # Fixed RNA Profiling (Singleplex)
    MFRP = "MFRP"  # Fixed RNA Profiling (Multiplex, Probe Barcode on R2)
    MFRP_R1 = "MFRP-R1"  # Fixed RNA Profiling (Multiplex, Probe Barcode on R1)
    MFRP_RNA = "MFRP-RNA"  # Fixed RNA Profiling (Multiplex, RNA, Probe Barcode on R2)
    MFRP_Ab = "MFRP-Ab"  # Fixed RNA Profiling (Multiplex, Antibody, Probe Barcode at R2:69)
    MFRP_Ab_R2pos50 = "MFRP-Ab-R2pos50"  # Fixed RNA Profiling (Multiplex, Antibody, Probe Barcode at R2:50)
    MFRP_RNA_R1 = "MFRP-RNA-R1"  # Fixed RNA Profiling (Multiplex, RNA, Probe Barcode on R1)
    MFRP_Ab_R1 = "MFRP-Ab-R1"  # Fixed RNA Profiling (Multiplex, Antibody, Probe Barcode on R1)
    ARC_v1 = "ARC-v1"  # for analyzing the Gene Expression portion of Multiome data. If Cell Ranger auto-detects ARC-v1 chemistry, an error is triggered.


five_prime_seq = typer.Typer(
    name="5'-scRNAseq script generator",
    help=(
        "Use information from a bcl-convert samplesheet to create scripts to process data from the 10x Genomics"
        " Single Cell Immune Profiling data using cellranger",
    ),
    add_completion=False
)


@five_prime_seq.command(name="version", context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def version_callback(value: Annotated[bool, typer.Option()] = True) -> None:  # FBT001
    """Prints the version of the package."""
    if value:
        rprint(f"[yellow]asapseq-script-gen[/] version: [bold blue]{__version__}[/]")
        raise typer.Exit()


def create_5_prime_multi_body(
    sample: str,
    multi_samplesheet: Path,
    job_interval: int = 2000,
    max_num_jobs: int = 8,
    mem_per_core: int = 32,
    job_template_path: Path = JOB_MANAGER_TEMPLATE_PATH,
) -> str:
    return (
        f"cellranger \\\n"
        f"\tmulti \\\n"
        f"\t\t--id={sample} \\\n"
        f"\t\t--csv={multi_samplesheet} \\\n"
        f"\t\t--jobinterval={job_interval} \\\n"
        f"\t\t--jobmode={job_template_path} \\\n"
        f"\t\t--mempercore={mem_per_core} \\\n"
        f"\t\t--maxjobs={max_num_jobs}"
    )


def create_gex_block(
    create_bam: bool = False,
    expected_cells: int = 20000,
    include_introns: bool = True,
    no_secondary_analysis: bool = True,
    chemistry: Literal["SC5PHT"] = "SC5PHT",  # This is only here for future expansion
    reference: Path = GEX_REFERENCE_PATH,
) -> str:
    return (
        "[gene-expression],,,\n"
        f"create-bam,{str(create_bam).lower()},,\n"
        f"expect-cells,{expected_cells},,\n"
        f"include-introns,{str(include_introns).lower()},,\n"
        f"no-secondary,{str(no_secondary_analysis).lower()},,\n"
        f"chemistry,{chemistry},,\n"
        f"reference,{reference},,\n"
    )


def create_vdj_block(
    reference: Path = VDJ_REFERENCE_PATH,
) -> str:
    return f"[vdj],,,\nreference,{reference},,\n"


def create_feature_block(
    reference: Path = FEAT_REFERENCE_PATH,
) -> str:
    return f"[feature],,,\nreference,{reference},,\n"


def create_lib_block(
    sample: str,
    fastq_path: Path,
    gex: bool = True,
    use_gex_lower: bool = True,
    gex_lanes: int | list[int] | None = None,
    bcr: bool = True,
    use_bcr_lower: bool = True,
    bcr_lanes: int | list[int] | None = None,
    tcr: bool = True,
    use_tcr_lower: bool = False,
    tcr_lanes: int | list[int] | None = None,
    feat: bool = True,
    use_feat: bool = True,
    feat_lanes: int | list[int] | None = None,
    feat_suffix: str = "prot",
) -> str:
    header = "[libraries],,,\nfastq_id,fastqs,lanes,feature_types\n"
    gex_row = (
        f"{sample}_{'gex' if use_gex_lower else 'GEX'},{fastq_path},{gex_lanes or ''},Gene Expression\n"
        if gex
        else ""
    )
    bcr_row = (
        f"{sample}_{'bcr' if use_bcr_lower else 'BCR'},{fastq_path},{bcr_lanes or ''},VDJ-B\n" if bcr else ""
    )
    tcr_row = (
        f"{sample}_{'tcr' if use_tcr_lower else 'TCR'},{fastq_path},{tcr_lanes or ''},VDJ-T\n" if tcr else ""
    )
    feat_row = (
        f"{sample}_{'feat' if use_feat else 'prot'},{fastq_path},{feat_lanes or ''},Antibody Capture\n"
        if feat
        else ""
    )

    return header + gex_row + bcr_row + tcr_row + feat_row


def create_multi_samplesheet(
    sample: str,
    fastq_path: Path,
    output: Path,
    create_bam: bool = False,
    expected_cells: int = 20000,
    include_introns: bool = True,
    no_secondary_analysis: bool = True,
    chemistry: str = "SC5PHT",
    gex_reference: Path = GEX_REFERENCE_PATH,
    vdj_reference: Path = VDJ_REFERENCE_PATH,
    feat_reference: Path = FEAT_REFERENCE_PATH,
    gex: bool = True,
    use_gex_lower: bool = True,
    bcr: bool = True,
    use_bcr_lower: bool = True,
    tcr: bool = True,
    use_tcr_lower: bool = True,
    feat: bool = True,
    use_feat: bool = True,
) -> Path:
    samplesheet_data = (
        create_gex_block(
            create_bam=create_bam,
            expected_cells=expected_cells,
            include_introns=include_introns,
            no_secondary_analysis=no_secondary_analysis,
            chemistry=chemistry,
            reference=gex_reference,
        )
        + ",,,\n"
        + create_vdj_block(
            reference=vdj_reference,
        )
        + ",,,\n"
        + create_feature_block(reference=feat_reference)
        + ",,,\n"
        + create_lib_block(
            sample=sample,
            fastq_path=fastq_path,
            gex=gex,
            use_gex_lower=use_gex_lower,
            bcr=bcr,
            use_bcr_lower=use_bcr_lower,
            tcr=tcr,
            use_tcr_lower=use_tcr_lower,
            feat=feat,
            use_feat=use_feat,
        )
    )

    with output.open("w") as s:
        s.writelines(samplesheet_data)

    return output


@five_prime_seq.command(name="create_scripts", no_args_is_help=True)
def create_five_prime_script(
    samplesheet: Annotated[
        Path,
        typer.Argument(help="Path to the bcl-convert samplesheet"),
    ],
    scripts_out_folder: Annotated[
        Path,
        typer.Argument(help="Path to where the generated scripts should be written"),
    ],
    fastq_path: Annotated[Path, typer.Argument(help="Path to where the FASTQs are stored")],
    gex_index_path: Annotated[
        Path,
        typer.Option("--gex_index", help="Path to the salmon index of the CITE-seq barcodes"),
    ] = GEX_REFERENCE_PATH,
    vdj_index_path: Annotated[
        Path,
        typer.Option("--vdj_index", help="Path to the salmon index of the CITE-seq barcodes"),
    ] = VDJ_REFERENCE_PATH,
    feat_index_ref: Annotated[
        Path,
        typer.Option("--feat_ref", help="Path to the salmon index of the CITE-seq barcodes"),
    ] = FEAT_REFERENCE_PATH,
    kit_chemistry: Annotated[
        Optional[Chemistry], typer.Option("--chem", help="10x kit chemistry. Currently, only 'SC5PHT' is supported")
    ] = Chemistry.SC5PHT,
    job_interval: Annotated[
        int,
        typer.Option("--interval"),
    ] = 2000,
    max_num_jobs: Annotated[
        int,
        typer.Option("--max_jobs"),
    ] = 8,
    mem_per_core: Annotated[
        int,
        typer.Option("--mem_per_core"),
    ] = 8,
    mem: Annotated[
        int,
        typer.Option("--mem", "-m", help="Amount of memory to use for each slurm job"),
    ] = 64,
    cpus: Annotated[
        int,
        typer.Option("--cpus", "-c", help="Amount of cpus to use for each slurm job"),
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
        ),
    ] = False,
):
    ss_df = read_samplesheet(samplesheet)

    if "Sample_Project" not in ss_df.columns:
        msg = "'Sample_Project' was not found in the samplesheet columns. Is your samplesheet in the correct format?"
        raise ColumnNotFoundError(msg)

    module_line = "module load cellranger/8.0.0\n" if load_cellranger_module else "\n"

    libraries = {
        i: {
            j: ss_df.filter(ss_df["Sample_Project"] == i)["Sample_ID"].str.ends_with(j).any()
            for j in ["GEX", "gex", "BCR", "bcr", "TCR", "tcr", "prot", "feat"]
        }
        for i in ss_df["Sample_Project"].unique()
    }

    for k in libraries:
        gex_libraries_present = libraries[k]["GEX"] or libraries[k]["gex"]
        bcr_libraries_present = libraries[k]["BCR"] or libraries[k]["bcr"]
        tcr_libraries_present = libraries[k]["TCR"] or libraries[k]["tcr"]
        feat_libraries_present = libraries[k]["prot"] or libraries[k]["feat"]

        part1 = create_script_header(sample=k, jobtype="multi_count", mem=mem, cpus=cpus)
        part2 = create_multi_samplesheet(
            sample=k,
            fastq_path=fastq_path.joinpath(k),
            output=scripts_out_folder.joinpath(f"{k}_multi_samplesheet.csv").absolute(),
            create_bam=False,
            expected_cells=2000,
            include_introns=True,
            no_secondary_analysis=True,
            chemistry=kit_chemistry,
            gex_reference=gex_index_path,
            vdj_reference=vdj_index_path,
            feat_reference=feat_index_ref,
            gex=gex_libraries_present,
            use_gex_lower=libraries[k]["gex"],
            bcr=bcr_libraries_present,
            use_bcr_lower=libraries[k]["bcr"],
            tcr=tcr_libraries_present,
            use_tcr_lower=libraries[k]["tcr"],
            feat=feat_libraries_present,
            use_feat=libraries[k]["feat"],
        )
        part3 = create_5_prime_multi_body(
            sample=k,
            multi_samplesheet=part2,
            job_interval=job_interval,
            max_num_jobs=max_num_jobs,
            mem_per_core=mem_per_core,
            job_template_path=job_template_path,
        )
        script_text = part1 + module_line + part3
        outfile = scripts_out_folder.joinpath(f"{k}_multi_count.job")
        with outfile.open("w") as sf:
            sf.writelines(script_text)
