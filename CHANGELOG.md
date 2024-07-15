# CHANGELOG

# [2.4.0] - 2024-07-12

## Changed

- `sc_script_gen.scrnaseq.create_scrnaseq_script` can now take additional options to alter what is put in the 
    multi_samplesheet.csv file
- Changed the default value for the `expect_cells` argument to `None`

# [2.3.0] - 2024-06-11

## Fixed

- `sc_script_gen.utils.read_samplesheet` can now process samplesheets for bcl2fastq

## Changed

- Renamed the `fiveprime` subcommand `scrnaseq` since it covered more than just the
    five prime assays

# [2.2.0] - 2024-06-17

## Added

- Progress bars! Feedback!

# [2.1.1] - 2024-06-17

## Fixed

- Gee, maybe it would be nice if the ASAPseq protein counts output from Salmon didn't overwrite every other Salmon job?
    Now outputs to a sample-specific directory

# [2.1.0] - 2024-06-17

## Added

- README.md

## Fixed

- changed the `libtype` argument to salmon alevin in `create_salmon_script_body` to the proper `libType`

# [2.0.3] - 2024-06-17

## Fixed

- changed the `libtype` argument to salmon alevin in `create_salmon_script_body` to the proper `libType`

# [2.0.2] - 2024-06-10

## Fixed

- scrnaseq multi_samplesheet is no longer created with a duplicate sample name (gods, I need to start writing 
    the unit tests first)

# [2.0.1] - 2024-06-10

## Changed

- `job_interval` now defaults to the proper 2000 ms, not 5
- `mem` now defaults to 64

## Fixed

- multi_samplesheet is defined using an absolute path now
- gex/vdj/feat_indices now actually use the constant defined at the top of the module and do not default to `None`

# [2.0.0] - 2024-06-10

## Added

- Now writes scripts for Single Cell Immune Profiling (i.e. 5' assays)

## Changed

- Renamed module
- Renamed `create_count_scripts` submodule to `asapseq`
- Moved some shared functions to `utils` submodule

# [1.1.1] - 2024-06-07

## Fix

- Hotfix for "unterminated string literal"

# [1.1.0] - 2024-06-07

## Changed

- reformated code calling subfunctions to generate script parts so as to make arguments explicit

## Fixed

- Fixed a spot where a couple or arguments to a subfunction were mixed up
- Changed a few double quotes within a double-quote-enclosed f-string to single quote so this works with 
    Python < 3.12

# [1.0.0] - 2024-06-07

## Added

- Everything

[2.4.0]: https://github.com/milescsmith/asap-script-gen/compare/2.3.0..2.4.0
[2.3.0]: https://github.com/milescsmith/asap-script-gen/compare/2.2.0..2.3.0
[2.2.0]: https://github.com/milescsmith/asap-script-gen/compare/2.1.1..2.2.0
[2.1.1]: https://github.com/milescsmith/asap-script-gen/compare/2.1.0..2.1.1
[2.1.0]: https://github.com/milescsmith/asap-script-gen/compare/2.0.3..2.1.0
[2.0.3]: https://github.com/milescsmith/asap-script-gen/compare/2.0.2..2.0.3
[2.0.2]: https://github.com/milescsmith/asap-script-gen/compare/2.0.1..2.0.2
[2.0.1]: https://github.com/milescsmith/asap-script-gen/compare/2.0.0..2.0.1
[2.0.0]: https://github.com/milescsmith/asap-script-gen/compare/1.1.1..2.0.0
[1.1.1]: https://github.com/milescsmith/asap-script-gen/compare/1.1.0..1.1.1
[1.1.0]: https://github.com/milescsmith/asap-script-gen/compare/1.0.0..1.1.0
[1.0.0]: https://github.com/milescsmith/asap-script-gen/releases/tag/1.0.0