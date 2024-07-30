from pathlib import Path

from sc_script_gen.scrnaseq import create_scrnaseq_script

create_scrnaseq_script(
    samplesheet=Path("/mnt/scratch/guth-aci/ALE06/metadata/20230907_set_1/20230907_ALE06_set_1_pool_1_samplesheet.csv"),
    scripts_out_folder=Path("/mnt/scratch/guth-aci/ALE06/scripts/20230907_set_1/counts"),
    fastq_path=Path("/mnt/scratch/guth-aci/ALE06/data/raw/fastqs/set_1_pool_1"),
)
