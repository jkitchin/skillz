# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial release of Claude Skills CLI tool
- Complete CLI for managing skills and commands
- 6 initial skills:
  - `python-ase`: Atomic Simulation Environment helper
  - `python-best-practices`: Python coding standards
  - `emacs-lisp`: Emacs Lisp development assistant
  - `fairchem`: Fair chemistry tools
  - `eln`: Electronic lab notebook management
  - `phd-qualifier`: PhD qualifier review
- Commands: install, uninstall, list, search, info, update, create
- Validation system for skills and commands
- Multi-platform support (Claude, Codex, Gemini)
- Configuration management
- Template system for creating new skills/commands
- Comprehensive test suite (44 tests)

## [0.1.0] - 2024-11-05

### Added
- Initial project structure
- Basic CLI framework with Click
- Configuration system with YAML
- Skill and command validators
- File operations utilities
- Rich terminal output
- MIT License
- Documentation (README.md, CONTRIBUTING.md)

[Unreleased]: https://github.com/jkitchin/skillz/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/jkitchin/skillz/releases/tag/v0.1.0
