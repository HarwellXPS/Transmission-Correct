# HarwellXPS — VAMAS Transmission-Function Correction Tool

A small, dependency-light tool for **transmission-function correction of XPS spectra** stored in the ISO-14976 **VAMAS** (`.vms`) format, as exported by **CasaXPS** and similar instruments.

It divides each spectrum by its encoded transmission / intensity-response function (`I → I / T`) and writes out a corrected VAMAS file plus a side-by-side comparison CSV. It works on a single file through an interactive viewer, or on an entire folder of files at once — both from a GUI and from the command line.

![Python](https://img.shields.io/badge/python-3.8%2B-blue)
![Status](https://img.shields.io/badge/status-stable-brightgreen)
![License](https://img.shields.io/badge/license-MIT-lightgrey)

---

## Purpose

This tool is aimed at two related needs:

1. **Permanently apply the data's own transmission function.** It bakes the instrument
   transmission / intensity-response function — exactly as encoded in the VAMAS file —
   into the spectra, writing a corrected `.vms` in which the intensities are already
   `I / T` and the transmission variable is reset to `1.0`. The correction is then
   "built in" and will not be applied a second time on re-import.
2. **Produce data sets for plotting the transmission-corrected data.** The comparison
   CSV holds the raw and corrected intensities side by side, ready for plotting in
   Excel, Origin, or pandas.

It performs transmission correction only — not background subtraction, peak fitting,
or quantification.

---

## Screenshots

> Replace the placeholders below with your own captures. Put the images in
> `docs/screenshots/` (create the folder) and the links will resolve on GitHub.

**Interactive viewer — raw vs. corrected**

<!-- ![Interactive viewer](docs/screenshots/viewer.png) -->
<p align="center"><em>docs/screenshots/viewer.png — add a capture of the main window</em></p>

**Bulk correction dialog**

<!-- ![Bulk correction](docs/screenshots/bulk.png) -->
<p align="center"><em>docs/screenshots/bulk.png — add a capture of the bulk folder dialog</em></p>

---

## Features

- **Transmission correction** — applies `I / T` per energy channel; channels with `T ≤ 0` are set to zero.
- **Interactive viewer** — load a `.vms` file and inspect raw vs. corrected spectra side-by-side or stacked, with a Binding-Energy / Kinetic-Energy axis toggle and multi-block selection.
- **Bulk processing** — correct every `.vms` file in a folder (optionally recursing into subfolders) via a threaded GUI dialog with a live progress bar and log, or headless from the CLI.
- **Two output formats per file**
  - a **corrected `.vms`** (counts replaced by `I / T`, transmission reset to `1.0`)
  - a **comparison `.csv`** containing both raw and corrected intensities for easy plotting in Excel, Origin, or pandas.
- **Selective export** — save *all* blocks (correcting only the selected ones) or save *only* the selected regions as a new, valid VAMAS file with its header block count updated.
- **Metadata-aware writer** — recomputes the per-variable min/max ordinate values when saving, so downstream readers (CasaXPS, etc.) autoscale the corrected data correctly.
- **Resilient batches** — a malformed file is logged and skipped rather than aborting the whole run.

---

## Requirements

| Component    | Needed for          | Notes                                  |
|--------------|---------------------|----------------------------------------|
| Python 3.8+  | everything          |                                        |
| `numpy`      | everything          | core numerics                          |
| `matplotlib` | GUI only            | interactive plots                      |
| `tkinter`    | GUI only            | ships with most Python installs        |

The headless bulk CLI needs only **Python + NumPy** — no display, no matplotlib, no Tk.

---

## Usage

### Interactive GUI (single file)

```bash
python xps_transmission_tool.py
```

1. **Load VAMAS File** — open a `.vms` file; detected blocks appear in the left panel.
2. Select one or more blocks (Ctrl/Shift-click, or **Select All**).
3. Toggle **Binding / Kinetic Energy** and the **Side-by-Side / Stacked** layout to inspect raw vs. corrected.
4. **Correct Selected & Save VMS** — write a corrected VAMAS file for the chosen blocks.
   - By default the output keeps **all** blocks and corrects only the selected ones; the rest are copied through unchanged. A confirmation dialog spells this out before saving.
   - Tick **Save only selected blocks** to instead export a VAMAS file containing **only** the selected regions (all corrected), with the file's "number of blocks" header field updated automatically. If the file's block boundaries can't be verified, the tool says why and offers to save all blocks instead.
5. **Export Comparison CSV** — write a CSV with raw and corrected intensities.
6. **Bulk Correct Folder…** — process a whole directory (see below) with progress and a log.

### Bulk correction (command line, headless)

```bash
python xps_transmission_tool.py --bulk INPUT_DIR [options]
```

| Option         | Description                                              |
|----------------|----------------------------------------------------------|
| `--out DIR`    | Output directory (default: `INPUT_DIR/corrected`).       |
| `--recurse`    | Also process `.vms` files in subdirectories.             |
| `--no-vms`     | Skip writing corrected `.vms` files.                     |
| `--no-csv`     | Skip writing comparison CSV files.                       |

**Examples**

```bash
# Correct every .vms in ./data, write results to ./data/corrected
python xps_transmission_tool.py --bulk ./data

# Recurse, send output elsewhere, CSV only
python xps_transmission_tool.py --bulk ./data --out ./results --recurse --no-vms
```

The CLI prints a per-file `OK`/`FAIL` log and returns a non-zero exit code if any file failed — convenient for scripting.

---

## Output files

For each input `sample.vms`:

| File                      | Contents                                                                 |
|---------------------------|--------------------------------------------------------------------------|
| `sample_corrected.vms`    | VAMAS file with counts replaced by `I / T` and transmission reset to `1.0`. Structure and all other blocks preserved. |
| `sample_comparison.csv`   | Long-format table for plotting/analysis.                                  |

The comparison CSV has one row per energy channel:

```csv
Region,Binding_Energy_eV,Kinetic_Energy_eV,Raw_Intensity_counts,Transmission,Corrected_Intensity
[1] Survey,1386.600,100.000,1000.0,0.500000,2.000000e+03
...
```

The `Region` column lets multiple blocks (which may have different energy axes) coexist in a single file — filter by region in Excel/pandas/Origin.

---

## How the correction works

Each VAMAS data point stores two corresponding variables in the CasaXPS convention:

1. **counts** — the measured intensity, and
2. **transmission** — the spectrometer's transmission / intensity-response function.

The tool computes the corrected intensity as

```
I_corrected = counts / transmission
```

and, in the written `.vms`, resets the transmission variable to `1.0` so the correction is not applied a second time on re-import. The per-variable minimum/maximum ordinate metadata is recomputed to match the corrected data.

---

## Configuration & assumptions

A couple of constants near the top of `xps_transmission_tool.py` control instrument-specific behaviour:

- **`N_CORR_VARS = 2`** — each point is assumed to carry two corresponding variables, *counts* then *transmission*. Change this (and the column mapping) if your exports differ.
- **`REGION_NAMES`** — only a *fallback*. Region names are normally read straight from the VAMAS experiment header and each block's identifier (e.g. `Ce 3d/2`), so this list is used only when the header can't be parsed.
- **`SOURCE_LABELS` / `DEFAULT_EXCITATION_EV`** — only a *fallback*. The photon energy for the KE ↔ BE conversion is normally read per block from the file's analysis-source field, so mixed-anode files (e.g. Al + Ag) convert correctly without configuration.

---

The code is organised as a GUI-free **core** (parser → correction → writers → bulk driver) with the Tkinter/matplotlib interface layered on top, so the parsing and correction logic can be imported and reused independently:

```python
from xps_transmission_tool import VamasParser, write_corrected_vms, write_comparison_csv

p = VamasParser("sample.vms")
indices = [b.index for b in p.blocks]
write_corrected_vms(p, indices, "sample_corrected.vms")
write_comparison_csv(p, indices, "sample_comparison.csv")
```

---

## Limitations & notes

- The VAMAS block parser uses pragmatic field-offset heuristics tuned to CasaXPS-style exports. If a file returns **zero blocks** or unexpected point counts, the field layout likely differs — open an issue with a sample file.
- Only the *counts* and *transmission* corresponding variables are used; additional variables are not currently handled.
- This tool performs transmission correction only. It does **not** do background subtraction, peak fitting, or quantification.

---

## Contributing

Issues and pull requests are welcome. If you hit a file that doesn't parse, attaching a (small, anonymised) sample `.vms` is the fastest way to get support added.

---

## License

Released under the MIT License. See [`LICENSE`](LICENSE) for details.

---

## Acknowledgements

Built for streamlining transmission-function correction of XPS VAMAS exports. VAMAS is the ISO-14976 surface-chemical-analysis data transfer format.
