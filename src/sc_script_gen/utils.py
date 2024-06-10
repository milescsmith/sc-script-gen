from pathlib import Path

import polars as pl

JOB_MANAGER_TEMPLATE_PATH = "/Volumes/guth_aci_informatics/software/slurm.template"


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
) -> str:
    """Generate a standard slurm job file header"""
    return f"#! /bin/bash -l\n\n#SBATCH -J {sample}{f'_{jobtype}' if jobtype else ''}\n#SBATCH -o {sample}{f'_{jobtype}' if jobtype else ''}.log\n#SBATCH --mail-user=None\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mem={mem}G\n#SBATCH --partition=serial\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task={cpus}\n"
