from pathlib import Path
from sc_script_gen.scrnaseq import create_five_prime_script

create_five_prime_script(
    samplesheet=Path("/mnt/scratch/guth-aci/ARA08/metadata/samplesheets/CGC-MS-2096_20240410_ara08_scRNAseq_set_2.csv"),
    scripts_out_folder=Path("/home/milo/workspace/sc_script_gen/tests/test_fiveprime"),
    fastq_path=Path("/mnt/scratch/guth-aci/ARA08/data/raw/fastqs/scrnaseq_set_2/"),
)