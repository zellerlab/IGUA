# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/zellerlab/IGUA/compare/v0.2.0...HEAD


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
