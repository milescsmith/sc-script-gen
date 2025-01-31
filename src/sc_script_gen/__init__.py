"""
asap-script-gen

Use the information in a a bcl-convert samplesheet
to create the slurm scripts necessary to process the data.

You are already demuxing the data, so you have the samplesheet -
why do all this again?
"""

from importlib.metadata import PackageNotFoundError, version

# from rich.console import Console

try:
    __version__ = version(__package__)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
