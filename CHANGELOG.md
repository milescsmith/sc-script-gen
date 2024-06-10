# CHANGELOG

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

[1.1.1]: https://github.com/milescsmith/asap-script-gen/compare/1.1.1..2.0.0
[1.1.1]: https://github.com/milescsmith/asap-script-gen/compare/1.1.0..1.1.1
[1.1.0]: https://github.com/milescsmith/asap-script-gen/compare/1.0.0..1.1.0
[1.0.0]: https://github.com/milescsmith/asap-script-gen/releases/tag/1.0.0