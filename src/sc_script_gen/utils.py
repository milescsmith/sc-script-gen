from pathlib import Path

import polars as pl
from loguru import logger

JOB_MANAGER_TEMPLATE_PATH = "/Volumes/guth_aci_informatics/software/slurm.template"


def read_samplesheet(samplesheet_file: Path) -> pl.DataFrame:
    """This should match either bcl-convert- or bcl2fastq-style samplesheets"""
    with samplesheet_file.open("r") as s:
        n = 0
        found = False
        for line in s.readlines():
            n += 1
            if "Data" in line:
                found = True
                break

    if not found:
        msg = "A line with 'Data' was not found in the samplesheet. Unable to process"
        logger.error(msg)
        raise ValueError(msg)
    return pl.read_csv(samplesheet_file, skip_rows=n)


def create_script_header(
    sample: str,
    jobtype: str | None = None,
    mem: int = 32,
    cpus: int = 1,
) -> str:
    """Generate a standard slurm job file header"""
    return f"#! /bin/bash -l\n\n#SBATCH -J {sample}{f'_{jobtype}' if jobtype else ''}\n#SBATCH -o {sample}{f'_{jobtype}' if jobtype else ''}.log\n#SBATCH --mail-user=None\n#SBATCH --mail-type=END,FAIL\n#SBATCH --mem={mem}G\n#SBATCH --partition=serial\n#SBATCH --nodes=1\n#SBATCH --cpus-per-task={cpus}\n"
