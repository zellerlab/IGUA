# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/zellerlab/IGUA/compare/v0.2.2...HEAD


## [v0.2.2] - 2026-08-12
[v0.2.2]: https://github.com/zellerlab/IGUA/compare/v0.2.1...v0.2.2

### Fixed
- Incorrect cluster length extraction for antiSMASH GenBank datasets ([#5](https://github.com/zellerlab/IGUA/issues/5)).


## [v0.2.1] - 2026-06-04
[v0.2.1]: https://github.com/zellerlab/IGUA/compare/v0.2.0...v0.2.1

### Fixed
- Missing import in `igua.cli.inputs` module ([#3](https://github.com/zellerlab/IGUA/issues/3)).

### Changed
- Update `pyo3` to `v0.28.3`.
- Relax `gb-io` dependency to support `v0.4`.


## [v0.2.0] - 2026-02-17
[v0.2.0]: https://github.com/zellerlab/IGUA/compare/v0.1.0...v0.2.0

### Added
- Public API to configure and run IGUA from Python code.
- Support for antiSMASH GenBank and Zip files.
- Option to disable weighing by protein length in distance matrix creation.
- Linear clustering method inspired by `mmseqs linclust` as an alternative to hierarchical clustering.

### Changed
- Rename flags in CLI to select the dataset type.
- Improve handling of temporary directory in CLI and pipeline code.


## [v0.1.0] - 2025-05-06
[v0.1.0]: https://github.com/zellerlab/IGUA/compare/6c5c7b...v0.1.0

Initial release.
