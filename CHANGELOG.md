# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and the project adheres to
[Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-06-13

### Added
- Per-block photon (excitation) energy is read directly from each block's source
  field, so mixed-anode files (e.g. Al + Ag) get correct binding-energy axes. The
  parser handles `MAP`/`MAPDP` experiment modes and Thermo as well as Kratos exports.
- Region names and block boundaries are read from the VAMAS experiment header and
  each block's identifier (e.g. `Ce 3d/2`), so files using arbitrary Kratos/CasaXPS
  region labels are named correctly and support selective export. `REGION_NAMES` is
  now only a fallback for files whose header can't be parsed.
- Interactive GUI viewer for VAMAS (`.vms`) XPS spectra: raw vs. transmission-corrected
  comparison, Binding/Kinetic energy axis toggle, side-by-side or stacked layouts, and
  multi-block selection.
- Transmission-function correction (`I / T`) producing a corrected `.vms` and a
  comparison `.csv` (raw and corrected intensities, long format with a `Region` column).
- Bulk correction of a whole directory, via a threaded GUI dialog with progress/log or
  the headless CLI (`--bulk`, with `--recurse`, `--out`, `--no-vms`, `--no-csv`).
- **Save only selected blocks** option: writes a valid VAMAS subset with the experiment
  header's "number of blocks" field updated and any trailer preserved. Guarded by
  block-boundary verification, with a graceful fallback if boundaries can't be confirmed.
- Confirmation warning when saving all blocks while only a subset is corrected.
- Per-variable min/max ordinate metadata recomputed on write so downstream readers
  autoscale corrected data correctly.
- Console entry point `harwell-xps`, packaging via `pyproject.toml`, a pytest suite, and
  GitHub Actions CI across Python 3.8-3.12.
