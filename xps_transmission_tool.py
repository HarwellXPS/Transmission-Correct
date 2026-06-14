#!/usr/bin/env python3
"""
HarwellXPS - VAMAS Transmission-Function Correction Tool
=========================================================

Reads ISO-14976 VAMAS (.vms) files exported by CasaXPS (and similar instruments),
divides each spectrum by its encoded transmission/intensity-response function, and
writes out:

  * a corrected VAMAS (.vms) file  (counts -> counts / T, transmission reset to 1.0)
  * a comparison CSV               (raw and corrected intensities side by side)

Intended use
------------
This tool is for people who want to **permanently bake the instrument transmission
function (as encoded in the data itself) into their spectra** -- producing files in
which the intensities are already transmission-corrected and the transmission
variable is reset to 1.0 so the correction is not applied again on re-import --
and/or to **generate data sets for plotting that show the transmission-corrected
data** (the comparison CSV holds raw and corrected intensities side by side). It
performs transmission correction only: it is not a background-subtraction,
peak-fitting, or quantification tool.

Two ways to run it:

  GUI (single file, interactive viewer):
      python xps_transmission_tool.py

  Headless bulk (correct every .vms in a folder):
      python xps_transmission_tool.py --bulk  INPUT_DIR  [--out OUTPUT_DIR]
                                      [--recurse] [--no-csv] [--no-vms]

Design notes / assumptions
--------------------------
* Each spectral block is assumed to carry TWO corresponding variables in the
  CasaXPS convention: variable 1 = counts, variable 2 = transmission function.
  This matches the original working parser. If your files differ, change
  ``N_CORR_VARS`` and the column mapping in ``VamasBlock``.
* Region/block names are read from the VAMAS experiment header and each block's
  identifier (e.g. ``Ce 3d/2``). If the header can't be parsed, the tool falls back
  to a known-name list (``REGION_NAMES``) with a ``Region N`` default.
* The photon (excitation) energy used for the KE<->BE conversion is read per block
  from its analysis-source field, so mixed-anode files (e.g. Al + Ag) convert
  correctly. ``SOURCE_LABELS`` / ``DEFAULT_EXCITATION_EV`` are only a fallback.
* On write-out, the per-variable min/max ordinate metadata is recomputed so that
  downstream readers (CasaXPS, etc.) autoscale the corrected data correctly --
  this was not handled by the original script.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from dataclasses import dataclass, field

import numpy as np

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

# Number of corresponding (ordinate) variables per data point.
# CasaXPS transmission-corrected exports interleave: counts, transmission.
N_CORR_VARS = 2

# Region identifiers recognised when naming blocks. Extend freely.
REGION_NAMES = {
    "Survey", "Valence",
    "Cu 2p", "Cu 2p/4", "C 1s", "O 1s", "Mo 3d", "S 2p",
    "N 1s", "F 1s", "Si 2p", "Ti 2p", "Fe 2p", "Au 4f",
    "Ag 3d", "Zn 2p", "Al 2p", "Na 1s", "Ca 2p", "Cl 2p",
}

# Source labels searched when recovering the photon (excitation) energy.
SOURCE_LABELS = {"Al mono", "Al", "Mg", "Al Ka", "Mg Ka"}
DEFAULT_EXCITATION_EV = 1486.6  # Al Kalpha, used if no source line is found.


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


# --------------------------------------------------------------------------- #
# Data model
# --------------------------------------------------------------------------- #

@dataclass
class VamasBlock:
    """A single parsed spectral region."""
    index: int                       # 0-based position in the file
    name: str                        # e.g. "[1] Survey"
    energy_type: str                 # "kinetic energy" | "binding energy"
    excitation_ev: float
    abscissa: np.ndarray             # energy axis as stored in the file
    ke: np.ndarray                   # kinetic energy axis
    be: np.ndarray                   # binding energy axis
    intensity: np.ndarray            # raw counts (corresponding variable 1)
    transmission: np.ndarray         # transmission function (corresponding variable 2)
    n_points: int
    data_start_line: int             # raw_lines index of first ordinate value
    minmax_start_line: int           # raw_lines index of first min/max metadata value
    block_start_line: int = -1       # raw_lines index of the block identifier (start)

    @property
    def data_end_line(self) -> int:
        """raw_lines index one past the block's last ordinate value."""
        return self.data_start_line + self.n_points * N_CORR_VARS

    def corrected(self) -> np.ndarray:
        """Transmission-corrected intensity: counts / T (0 where T <= 0)."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(self.transmission > 0,
                            self.intensity / self.transmission, 0.0)


# --------------------------------------------------------------------------- #
# Parser
# --------------------------------------------------------------------------- #

class VamasParser:
    """Dynamic state-machine parser for CasaXPS-style VAMAS files.

    Keeps the index arithmetic of the original tool (which is known to work on
    the target instrument exports) but isolates it behind a clean interface so
    it can be reused by both the GUI and the bulk pipeline.
    """

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.raw_lines: list[str] = []
        self.blocks: list[VamasBlock] = []
        self.header_end: int | None = None   # raw_lines index where block 0 begins
        self.count_line: int = -1            # raw_lines index of "number of blocks"
        self._exp_mode: str | None = None    # experiment mode (NORM, MAP, ...)
        self._n_exp_var: int = 0             # number of experimental variables
        self._parse()

    # -- public helpers ----------------------------------------------------- #

    def __len__(self) -> int:
        return len(self.blocks)

    def subset_info(self) -> dict:
        """Check whether a VAMAS file containing only a subset of blocks can be
        safely written.

        This requires (a) locating the experiment header's "number of blocks"
        field, which sits on the line immediately before the first block, and
        (b) confirming its value matches the number of parsed blocks. Returns
        ``{"ok": bool, "count_line": int, "reason": str}``.
        """
        if not self.blocks:
            return {"ok": False, "count_line": -1, "reason": "no blocks parsed"}

        start0 = self.blocks[0].block_start_line
        if start0 <= 0:
            return {"ok": False, "count_line": -1,
                    "reason": "could not locate the experiment header before the "
                              "first block (its region name may be unrecognised — "
                              "add it to REGION_NAMES)"}

        count_line = self.count_line if self.count_line >= 0 else start0 - 1
        try:
            declared = int(self.raw_lines[count_line])
        except (ValueError, IndexError):
            return {"ok": False, "count_line": -1,
                    "reason": "could not read the 'number of blocks' header field"}

        if declared != len(self.blocks):
            return {"ok": False, "count_line": -1,
                    "reason": f"header block count ({declared}) does not match the "
                              f"{len(self.blocks)} parsed blocks"}

        return {"ok": True, "count_line": count_line, "reason": ""}

    # -- internals ---------------------------------------------------------- #

    def _parse(self) -> None:
        with open(self.filepath, "r", errors="ignore") as fh:
            self.raw_lines = [line.strip() for line in fh.readlines()]

        lines = self.raw_lines
        n = len(lines)
        i = 0
        block_counter = 1

        while i < n:
            if lines[i].lower() in ("kinetic energy", "binding energy"):
                try:
                    block, next_i = self._parse_block(i, block_counter)
                    self.blocks.append(block)
                    block_counter += 1
                    i = next_i
                except Exception:
                    # Malformed / unexpected block: skip a line and keep going
                    # so one bad region never kills the whole file.
                    i += 1
            else:
                i += 1

        # Recover real region names and block boundaries from the experiment
        # header where possible (falls back to the REGION_NAMES heuristic).
        self._assign_block_identities()

    def _locate_first_block(self) -> tuple[int | None, int, int]:
        """Parse the experiment header to find where the first block starts.

        Returns ``(block0_start, count_line, declared_n_blocks)`` or
        ``(None, -1, 0)`` if the header cannot be parsed. ``count_line`` is the
        raw_lines index of the "number of blocks" field (the last header line).
        """
        lines = self.raw_lines
        try:
            j = 5                                   # 5 fixed identifier lines
            n_comment = int(lines[j]); j += 1 + n_comment
            mode = lines[j]; j += 1                 # experiment mode
            j += 1                                  # scan mode
            if mode in ("MAP", "MAPDP", "NORM", "SDP"):
                j += 1                              # number of spectral regions
            if mode in ("MAP", "MAPDP"):
                j += 3                              # analysis positions, nx, ny
            n_exp_var = int(lines[j]); j += 1 + 2 * n_exp_var
            n_incl_excl = int(lines[j]); j += 1 + abs(n_incl_excl)
            n_manual = int(lines[j]); j += 1 + n_manual
            n_future_exp = int(lines[j]); j += 1 + n_future_exp
            int(lines[j]); j += 1                   # number of future-upgrade block entries
            count_line = j
            declared = int(lines[count_line]); j += 1
            block0_start = j
            if not (0 < block0_start <= len(lines)):
                return None, -1, 0
            self._exp_mode = mode
            self._n_exp_var = n_exp_var
            return block0_start, count_line, declared
        except (ValueError, IndexError):
            return None, -1, 0

    def _assign_block_identities(self) -> None:
        """Override block names and start lines using the experiment header and
        VAMAS block contiguity, so real region identifiers (e.g. ``Ce 3d/2``)
        are used instead of the heuristic ``Region N`` fallback."""
        if not self.blocks:
            return

        block0_start, count_line, declared = self._locate_first_block()
        if block0_start is None or declared != len(self.blocks):
            return  # keep the REGION_NAMES fallback already assigned

        # Block k begins where block k-1's ordinate data ended (contiguous).
        starts = [block0_start] + [self.blocks[k - 1].data_end_line
                                   for k in range(1, len(self.blocks))]

        idents = []
        for s in starts:
            if not (0 <= s < len(self.raw_lines)):
                return
            ident = self.raw_lines[s]
            if not ident or _is_number(ident):
                return  # boundaries don't line up -> keep fallback
            idents.append(ident)

        self.header_end = block0_start
        self.count_line = count_line
        for k, blk in enumerate(self.blocks):
            blk.block_start_line = starts[k]
            blk.name = f"[{k + 1}] {idents[k]}"
            # Read the real photon energy from this block's source field so any
            # anode (Al, Mg, Ag, ...) gives a correct binding-energy axis.
            exc = self._read_block_excitation(starts[k])
            if exc is not None and exc > 0:
                blk.excitation_ev = exc
                if "kinetic" in blk.energy_type:
                    blk.ke = blk.abscissa
                    blk.be = exc - blk.abscissa
                else:
                    blk.be = blk.abscissa
                    blk.ke = exc - blk.abscissa

    def _read_block_excitation(self, block_start: int) -> float | None:
        """Read the analysis-source characteristic (photon) energy directly from
        a block header. Robust to any anode because it parses the structure
        rather than matching a fixed list of source labels."""
        lines = self.raw_lines
        try:
            idx = block_start + 9            # id, sample, date (6), GMT hours
            n_comment = int(lines[idx]); idx += 1 + n_comment
            idx += 1                         # technique
            if self._exp_mode in ("MAP", "MAPDP"):
                idx += 2                     # analysis x, y coordinates
            idx += self._n_exp_var           # experimental-variable values
            idx += 1                         # analysis source label
            return float(lines[idx])         # analysis source characteristic energy
        except (ValueError, IndexError):
            return None

    def _parse_block(self, i: int, block_counter: int) -> tuple[VamasBlock, int]:
        lines = self.raw_lines
        n = len(lines)

        energy_type = lines[i].lower()              # abscissa label
        start_energy = float(lines[i + 2])          # abscissa start value
        step_size = float(lines[i + 3])             # abscissa increment

        # Region name: backward trace for a known identifier.
        name, name_line = self._find_region_name(i)
        if not name:
            name = f"Region {block_counter}"
            name_line = -1

        # Advance to the signal-mode marker ("pulse counting").
        scan = i
        while scan < n and lines[scan] != "pulse counting":
            scan += 1
        if scan >= n:
            raise ValueError("signal-mode marker not found")

        idx = scan + 7                              # -> number of additional params
        num_supp = int(lines[idx])
        idx += 1 + num_supp * 3                     # skip additional (label,unit,value)*N

        total_ordinates = int(lines[idx])
        idx += 1

        minmax_start = idx                          # per-var (min,max) pairs
        idx += 2 * N_CORR_VARS

        data_start = idx
        n_points = total_ordinates // N_CORR_VARS

        # Read interleaved ordinate columns: [v0, v1, v0, v1, ...]
        columns = [np.empty(n_points) for _ in range(N_CORR_VARS)]
        for p in range(n_points):
            for c in range(N_CORR_VARS):
                columns[c][p] = float(lines[idx])
                idx += 1

        intensity = columns[0]
        transmission = columns[1] if N_CORR_VARS > 1 else np.ones(n_points)

        # Energy axes.
        abscissa = start_energy + step_size * np.arange(n_points)
        excitation = self._find_excitation(i)
        if "kinetic" in energy_type:
            ke = abscissa
            be = excitation - ke
        else:                                       # binding-energy abscissa
            be = abscissa
            ke = excitation - be

        block = VamasBlock(
            index=block_counter - 1,
            name=f"[{block_counter}] {name}",
            energy_type=energy_type,
            excitation_ev=excitation,
            abscissa=abscissa,
            ke=ke,
            be=be,
            intensity=intensity,
            transmission=transmission,
            n_points=n_points,
            data_start_line=data_start,
            minmax_start_line=minmax_start,
            block_start_line=name_line,
        )
        return block, idx

    def _find_region_name(self, i: int) -> tuple[str, int]:
        lines = self.raw_lines
        lo = max(0, i - 100)
        for b in range(i, lo - 1, -1):
            if lines[b] in REGION_NAMES:
                return lines[b], b
        return "", -1

    def _find_excitation(self, i: int) -> float:
        lines = self.raw_lines
        n = len(lines)
        lo = max(0, i - 50)
        for b in range(i, lo - 1, -1):
            if lines[b] in SOURCE_LABELS and b + 1 < n:
                try:
                    return float(lines[b + 1])
                except ValueError:
                    continue
        return DEFAULT_EXCITATION_EV


# --------------------------------------------------------------------------- #
# Output writers
# --------------------------------------------------------------------------- #

def _apply_correction_in_place(lines: list[str], blk: VamasBlock) -> None:
    """Overwrite one block's ordinate data in ``lines`` with corrected values.

    Counts become ``I / T``, transmission is reset to 1.0, and the per-variable
    min/max ordinate metadata is recomputed so downstream readers autoscale the
    corrected trace correctly.
    """
    corr = blk.corrected()
    ds = blk.data_start_line
    for p in range(blk.n_points):
        lines[ds + p * N_CORR_VARS] = f"{corr[p]:.4f}"          # counts -> corrected
        if N_CORR_VARS > 1:
            lines[ds + p * N_CORR_VARS + 1] = "1.000000"        # transmission -> 1.0

    mm = blk.minmax_start_line
    if corr.size:
        lines[mm] = f"{float(np.min(corr)):.4f}"                # var0 min
        lines[mm + 1] = f"{float(np.max(corr)):.4f}"           # var0 max
    if N_CORR_VARS > 1:
        lines[mm + 2] = "1.000000"                             # var1 (T) min
        lines[mm + 3] = "1.000000"                             # var1 (T) max


def write_corrected_vms(parser: VamasParser, block_indices, save_path: str) -> int:
    """Write a VAMAS file with the chosen blocks transmission-corrected.

    ALL blocks are kept; only those in ``block_indices`` are corrected, the rest
    are copied verbatim. Returns the number of blocks corrected.
    """
    out = list(parser.raw_lines)
    wanted = set(block_indices)
    corrected_count = 0
    for blk in parser.blocks:
        if blk.index in wanted:
            _apply_correction_in_place(out, blk)
            corrected_count += 1

    with open(save_path, "w", encoding="ascii") as fh:
        fh.write("\n".join(out) + "\n")
    return corrected_count


def write_corrected_vms_subset(parser: VamasParser, block_indices, save_path: str) -> int:
    """Write a VAMAS file containing ONLY the selected blocks (all corrected).

    The experiment header is preserved with its "number of blocks" field updated,
    and any trailing content (e.g. an "end of experiment" terminator) is kept.
    Raises ValueError if the file's block boundaries cannot be verified.
    """
    info = parser.subset_info()
    if not info["ok"]:
        raise ValueError("cannot export only selected blocks — " + info["reason"])

    wanted = set(block_indices)
    selected = [b for b in parser.blocks if b.index in wanted]
    if not selected:
        raise ValueError("no blocks selected")

    # Corrected copy of every selected block.
    corrected = list(parser.raw_lines)
    for blk in selected:
        _apply_correction_in_place(corrected, blk)

    # Update the header's "number of blocks" field.
    corrected[info["count_line"]] = str(len(selected))

    start0 = parser.blocks[0].block_start_line
    header = corrected[:start0]

    # VAMAS blocks are contiguous: block k starts where block k-1's data ended.
    starts = [start0] + [parser.blocks[k - 1].data_end_line
                         for k in range(1, len(parser.blocks))]

    body: list[str] = []
    for k, blk in enumerate(parser.blocks):
        if blk.index in wanted:
            body.extend(corrected[starts[k]:blk.data_end_line])

    trailer = corrected[parser.blocks[-1].data_end_line:]

    out = header + body + trailer
    with open(save_path, "w", encoding="ascii") as fh:
        fh.write("\n".join(out) + "\n")
    return len(selected)


def write_comparison_csv(parser: VamasParser, block_indices, save_path: str) -> None:
    """Write a tidy (long-format) CSV with raw and corrected intensities.

    One row per energy channel, with a Region column so multiple blocks (which
    may have different energy axes) coexist cleanly in a single file.
    """
    wanted = set(block_indices)
    with open(save_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "Region",
            "Binding_Energy_eV",
            "Kinetic_Energy_eV",
            "Raw_Intensity_counts",
            "Transmission",
            "Corrected_Intensity",
        ])
        for blk in parser.blocks:
            if blk.index not in wanted:
                continue
            corr = blk.corrected()
            for p in range(blk.n_points):
                w.writerow([
                    blk.name,
                    f"{blk.be[p]:.3f}",
                    f"{blk.ke[p]:.3f}",
                    f"{blk.intensity[p]:.1f}",
                    f"{blk.transmission[p]:.6f}",
                    f"{corr[p]:.6e}",
                ])


# --------------------------------------------------------------------------- #
# Bulk pipeline (shared by CLI and GUI)
# --------------------------------------------------------------------------- #

def find_vms_files(input_dir: str, recurse: bool = False) -> list[str]:
    matches = []
    if recurse:
        for root, _dirs, files in os.walk(input_dir):
            for f in files:
                if f.lower().endswith(".vms"):
                    matches.append(os.path.join(root, f))
    else:
        for f in sorted(os.listdir(input_dir)):
            if f.lower().endswith(".vms"):
                matches.append(os.path.join(input_dir, f))
    return sorted(matches)


def process_file(src_path: str, out_dir: str,
                 write_vms: bool = True, write_csv: bool = True) -> dict:
    """Correct every block in one file. Returns a per-file result dict."""
    parser = VamasParser(src_path)
    base = os.path.splitext(os.path.basename(src_path))[0]
    result = {"file": src_path, "blocks": len(parser), "vms": None, "csv": None}

    if len(parser) == 0:
        raise ValueError("no valid VAMAS blocks found")

    all_indices = [b.index for b in parser.blocks]

    if write_vms:
        vms_path = os.path.join(out_dir, f"{base}_corrected.vms")
        write_corrected_vms(parser, all_indices, vms_path)
        result["vms"] = vms_path

    if write_csv:
        csv_path = os.path.join(out_dir, f"{base}_comparison.csv")
        write_comparison_csv(parser, all_indices, csv_path)
        result["csv"] = csv_path

    return result


def bulk_process_directory(input_dir: str, out_dir: str,
                           recurse: bool = False,
                           write_vms: bool = True,
                           write_csv: bool = True,
                           log_cb=None, progress_cb=None) -> dict:
    """Correct every .vms file in a directory.

    ``log_cb(str)`` and ``progress_cb(done, total)`` are optional callbacks so
    the GUI can show live progress; the CLI passes simple print-based versions.
    """
    def log(msg):
        if log_cb:
            log_cb(msg)

    os.makedirs(out_dir, exist_ok=True)
    files = find_vms_files(input_dir, recurse=recurse)
    total = len(files)
    summary = {"total": total, "ok": 0, "failed": 0, "results": [], "errors": []}

    if total == 0:
        log("No .vms files found.")
        if progress_cb:
            progress_cb(0, 0)
        return summary

    log(f"Found {total} .vms file(s).")
    for n, src in enumerate(files, start=1):
        rel = os.path.relpath(src, input_dir)
        try:
            res = process_file(src, out_dir, write_vms=write_vms, write_csv=write_csv)
            summary["ok"] += 1
            summary["results"].append(res)
            log(f"[{n}/{total}] OK   {rel}  ({res['blocks']} blocks)")
        except Exception as exc:  # noqa: BLE001 - report and continue the batch
            summary["failed"] += 1
            summary["errors"].append((src, str(exc)))
            log(f"[{n}/{total}] FAIL {rel}  -> {exc}")
        if progress_cb:
            progress_cb(n, total)

    log(f"Done. {summary['ok']} succeeded, {summary['failed']} failed.")
    return summary


# --------------------------------------------------------------------------- #
# GUI (imported lazily so headless/CLI use needs no display or matplotlib)
# --------------------------------------------------------------------------- #

def launch_gui() -> None:
    import queue
    import threading
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    from matplotlib.backends.backend_tkagg import (
        FigureCanvasTkAgg, NavigationToolbar2Tk)
    from matplotlib.figure import Figure

    class XPSViewerApp:
        def __init__(self, root):
            self.root = root
            self.root.title("HarwellXPS - VAMAS Transmission Correction")
            self.root.geometry("1300x850")
            self.parser: VamasParser | None = None
            self._build_ui()

        # -- layout --------------------------------------------------------- #
        def _build_ui(self):
            top = ttk.Frame(self.root, padding=10)
            top.pack(side=tk.TOP, fill=tk.X)

            ttk.Button(top, text="Load VAMAS File (.vms)",
                       command=self.load_file).pack(side=tk.LEFT, padx=5)
            ttk.Button(top, text="Bulk Correct Folder…",
                       command=self.open_bulk_dialog,
                       style="Accent.TButton").pack(side=tk.LEFT, padx=5)

            self.axis_var = tk.StringVar(value="BE")
            ttk.Radiobutton(top, text="Binding Energy (eV)", variable=self.axis_var,
                            value="BE", command=self.update_plot).pack(side=tk.LEFT, padx=10)
            ttk.Radiobutton(top, text="Kinetic Energy (eV)", variable=self.axis_var,
                            value="KE", command=self.update_plot).pack(side=tk.LEFT, padx=10)

            ttk.Label(top, text="Layout:").pack(side=tk.LEFT, padx=(15, 5))
            self.layout_var = tk.StringVar(value="Side-by-Side")
            layout_menu = ttk.Combobox(top, textvariable=self.layout_var,
                                       values=["Side-by-Side", "Stacked Top/Bottom"],
                                       state="readonly", width=18)
            layout_menu.pack(side=tk.LEFT, padx=5)
            layout_menu.bind("<<ComboboxSelected>>", lambda e: self.update_plot())

            ttk.Button(top, text="Correct Selected & Save VMS",
                       command=self.correct_selected_and_save_vms).pack(side=tk.RIGHT, padx=5)
            self.save_only_selected = tk.BooleanVar(value=False)
            ttk.Checkbutton(top, text="Save only selected blocks",
                            variable=self.save_only_selected).pack(side=tk.RIGHT, padx=5)
            ttk.Button(top, text="Export Comparison CSV",
                       command=self.export_current_data).pack(side=tk.RIGHT, padx=5)

            paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
            paned.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

            left = ttk.LabelFrame(paned, text="VAMAS Blocks (Ctrl/Shift-click for multi-select)",
                                  padding=5)
            scroll = ttk.Scrollbar(left, orient=tk.VERTICAL)
            self.block_listbox = tk.Listbox(left, yscrollcommand=scroll.set,
                                            selectmode=tk.EXTENDED, width=32)
            scroll.config(command=self.block_listbox.yview)
            scroll.pack(side=tk.RIGHT, fill=tk.Y)
            self.block_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            self.block_listbox.bind("<<ListboxSelect>>", lambda e: self.update_plot())
            paned.add(left, weight=1)

            btns = ttk.Frame(left)
            btns.pack(side=tk.BOTTOM, fill=tk.X, pady=(5, 0))
            ttk.Button(btns, text="Select All",
                       command=self._select_all).pack(fill=tk.X)

            right = ttk.Frame(paned)
            paned.add(right, weight=4)
            self.fig = Figure(figsize=(8, 6))
            self.canvas = FigureCanvasTkAgg(self.fig, master=right)
            self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            NavigationToolbar2Tk(self.canvas, right).update()

            self.status = tk.StringVar(value="Load a .vms file, or use Bulk Correct Folder.")
            ttk.Label(self.root, textvariable=self.status, relief=tk.SUNKEN,
                      anchor=tk.W, padding=4).pack(side=tk.BOTTOM, fill=tk.X)

        # -- helpers -------------------------------------------------------- #
        def _select_all(self):
            self.block_listbox.selection_set(0, tk.END)
            self.update_plot()

        def _selected_block_indices(self):
            return [self.parser.blocks[i].index
                    for i in self.block_listbox.curselection()]

        # -- actions -------------------------------------------------------- #
        def load_file(self):
            path = filedialog.askopenfilename(
                filetypes=[("VAMAS Format", "*.vms"), ("All Files", "*.*")])
            if not path:
                return
            self.parser = VamasParser(path)
            self.block_listbox.delete(0, tk.END)
            for blk in self.parser.blocks:
                self.block_listbox.insert(tk.END, blk.name)
            if self.parser.blocks:
                self.block_listbox.selection_set(0)
                self.block_listbox.focus_set()
                self.status.set(f"Loaded {len(self.parser)} block(s) from "
                                f"{os.path.basename(path)}")
                self.update_plot()
            else:
                self.status.set("No valid VAMAS blocks found.")
                messagebox.showwarning("Parse Warning",
                                       "No valid ISO-VAMAS records matched.")

        def update_plot(self):
            if not self.parser:
                return
            sel = self.block_listbox.curselection()
            if not sel:
                return
            blk = self.parser.blocks[sel[0]]
            use_be = self.axis_var.get() == "BE"
            x = blk.be if use_be else blk.ke
            x_label = "Binding Energy (eV)" if use_be else "Kinetic Energy (eV)"
            raw_y = blk.intensity
            corr_y = blk.corrected()

            self.fig.clf()
            if self.layout_var.get() == "Side-by-Side":
                ax_raw = self.fig.add_subplot(1, 2, 1)
                ax_corr = self.fig.add_subplot(1, 2, 2)
            else:
                ax_raw = self.fig.add_subplot(2, 1, 1)
                ax_corr = self.fig.add_subplot(2, 1, 2)

            ax_raw.plot(x, raw_y, color="darkblue", linewidth=1.3)
            ax_raw.set_title("Uncorrected Spectral Profile")
            ax_raw.set_ylabel("Raw Intensity (Counts)")
            ax_corr.plot(x, corr_y, color="crimson", linewidth=1.3)
            ax_corr.set_title("Transmission Corrected (I / T)")
            ax_corr.set_ylabel("Corrected Intensity (A.U.)")
            for ax in (ax_raw, ax_corr):
                ax.set_xlabel(x_label)
                ax.grid(True, linestyle=":", alpha=0.6)
                if use_be:
                    ax.set_xlim(max(x), min(x))
                else:
                    ax.set_xlim(min(x), max(x))

            count = f"({len(sel)} blocks selected)" if len(sel) > 1 else ""
            self.fig.suptitle(f"XPS Panel Analysis: {blk.name} {count}",
                              fontsize=11, fontweight="bold")
            self.fig.tight_layout()
            self.canvas.draw()

        def correct_selected_and_save_vms(self):
            if not self.parser or not self.parser.blocks:
                return
            indices = self._selected_block_indices()
            if not indices:
                messagebox.showwarning("Selection Required",
                                       "Select one or more blocks to correct.")
                return

            total = len(self.parser.blocks)
            n_sel = len(indices)
            only_selected = self.save_only_selected.get()

            if only_selected:
                info = self.parser.subset_info()
                if not info["ok"]:
                    fall_back = messagebox.askyesno(
                        "Cannot export only selected",
                        "The block boundaries in this file could not be verified, "
                        "so saving only the selected block(s) isn't possible.\n\n"
                        f"Reason: {info['reason']}.\n\n"
                        "Save ALL blocks instead (only the selected ones corrected)?")
                    if not fall_back:
                        return
                    only_selected = False

            if not only_selected:
                # Explicit warning: the file keeps every block, only some corrected.
                proceed = messagebox.askokcancel(
                    "Confirm save — all blocks included",
                    f"The saved file will contain ALL {total} block(s), but only the "
                    f"{n_sel} selected block(s) will be transmission-corrected. "
                    f"The other {total - n_sel} block(s) are copied through unchanged.\n\n"
                    "Tick \u201cSave only selected blocks\u201d if you want a file "
                    "containing just the corrected regions.\n\nContinue?")
                if not proceed:
                    return

            base = os.path.splitext(os.path.basename(self.parser.filepath))[0]
            suffix = "_selected_corrected.vms" if only_selected else "_corrected.vms"
            path = filedialog.asksaveasfilename(
                defaultextension=".vms",
                filetypes=[("VAMAS Data File", "*.vms"), ("All Files", "*.*")],
                initialfile=f"{base}{suffix}",
                title=f"Export {'only ' if only_selected else ''}"
                      f"{n_sel} selected block(s)")
            if not path:
                return
            try:
                if only_selected:
                    n = write_corrected_vms_subset(self.parser, indices, path)
                    summary = (f"Saved {n} selected block(s), all transmission-corrected.")
                else:
                    n = write_corrected_vms(self.parser, indices, path)
                    summary = (f"Saved all {total} block(s); {n} transmission-corrected, "
                               f"{total - n} unchanged.")
                self.status.set(f"{summary}  ->  {os.path.basename(path)}")
                messagebox.showinfo("Correction Complete",
                                    f"{summary}\n\nSaved to:\n{path}")
            except Exception as exc:  # noqa: BLE001
                messagebox.showerror("Export Error", str(exc))

        def export_current_data(self):
            if not self.parser:
                return
            indices = self._selected_block_indices()
            if not indices:
                messagebox.showwarning("Export Warning",
                                       "Select one or more blocks to export.")
                return
            base = os.path.splitext(os.path.basename(self.parser.filepath))[0]
            path = filedialog.asksaveasfilename(
                defaultextension=".csv",
                filetypes=[("CSV", "*.csv"), ("All Files", "*.*")],
                initialfile=f"{base}_comparison.csv",
                title="Export comparison CSV")
            if not path:
                return
            try:
                write_comparison_csv(self.parser, indices, path)
                self.status.set(f"Exported CSV -> {os.path.basename(path)}")
                messagebox.showinfo("Export Complete", f"Saved to:\n{path}")
            except Exception as exc:  # noqa: BLE001
                messagebox.showerror("Export Error", str(exc))

        # -- bulk dialog ---------------------------------------------------- #
        def open_bulk_dialog(self):
            dlg = tk.Toplevel(self.root)
            dlg.title("Bulk Correct Folder")
            dlg.geometry("640x560")
            dlg.transient(self.root)

            state = {"in": tk.StringVar(), "out": tk.StringVar()}
            recurse = tk.BooleanVar(value=False)
            want_vms = tk.BooleanVar(value=True)
            want_csv = tk.BooleanVar(value=True)

            def pick_in():
                d = filedialog.askdirectory(title="Input folder of .vms files")
                if d:
                    state["in"].set(d)
                    if not state["out"].get():
                        state["out"].set(os.path.join(d, "corrected"))

            def pick_out():
                d = filedialog.askdirectory(title="Output folder")
                if d:
                    state["out"].set(d)

            frm = ttk.Frame(dlg, padding=12)
            frm.pack(fill=tk.BOTH, expand=True)

            ttk.Label(frm, text="Input folder:").grid(row=0, column=0, sticky=tk.W)
            ttk.Entry(frm, textvariable=state["in"], width=52).grid(row=0, column=1, padx=5)
            ttk.Button(frm, text="Browse…", command=pick_in).grid(row=0, column=2)

            ttk.Label(frm, text="Output folder:").grid(row=1, column=0, sticky=tk.W, pady=(8, 0))
            ttk.Entry(frm, textvariable=state["out"], width=52).grid(row=1, column=1, padx=5, pady=(8, 0))
            ttk.Button(frm, text="Browse…", command=pick_out).grid(row=1, column=2, pady=(8, 0))

            opts = ttk.Frame(frm)
            opts.grid(row=2, column=0, columnspan=3, sticky=tk.W, pady=10)
            ttk.Checkbutton(opts, text="Search subfolders", variable=recurse).pack(side=tk.LEFT, padx=5)
            ttk.Checkbutton(opts, text="Write corrected VMS", variable=want_vms).pack(side=tk.LEFT, padx=5)
            ttk.Checkbutton(opts, text="Write comparison CSV", variable=want_csv).pack(side=tk.LEFT, padx=5)

            progress = ttk.Progressbar(frm, length=560, mode="determinate")
            progress.grid(row=3, column=0, columnspan=3, pady=(4, 6))

            log_box = tk.Text(frm, height=16, width=74, state=tk.DISABLED, wrap=tk.NONE)
            log_box.grid(row=4, column=0, columnspan=3)

            events: "queue.Queue" = queue.Queue()

            def append_log(msg):
                log_box.config(state=tk.NORMAL)
                log_box.insert(tk.END, msg + "\n")
                log_box.see(tk.END)
                log_box.config(state=tk.DISABLED)

            def worker(in_dir, out_dir, rec, vms, csv_):
                bulk_process_directory(
                    in_dir, out_dir, recurse=rec, write_vms=vms, write_csv=csv_,
                    log_cb=lambda m: events.put(("log", m)),
                    progress_cb=lambda d, t: events.put(("prog", d, t)))
                events.put(("done", None))

            def poll():
                try:
                    while True:
                        ev = events.get_nowait()
                        if ev[0] == "log":
                            append_log(ev[1])
                        elif ev[0] == "prog":
                            done, total = ev[1], ev[2]
                            progress["maximum"] = max(total, 1)
                            progress["value"] = done
                        elif ev[0] == "done":
                            run_btn.config(state=tk.NORMAL)
                            self.status.set("Bulk correction finished.")
                            return
                except queue.Empty:
                    pass
                dlg.after(100, poll)

            def run():
                in_dir = state["in"].get().strip()
                out_dir = state["out"].get().strip()
                if not in_dir or not os.path.isdir(in_dir):
                    messagebox.showwarning("Input Required", "Choose a valid input folder.")
                    return
                if not out_dir:
                    messagebox.showwarning("Output Required", "Choose an output folder.")
                    return
                if not (want_vms.get() or want_csv.get()):
                    messagebox.showwarning("Nothing to do",
                                           "Enable VMS and/or CSV output.")
                    return
                run_btn.config(state=tk.DISABLED)
                append_log(f"Input : {in_dir}")
                append_log(f"Output: {out_dir}")
                threading.Thread(
                    target=worker,
                    args=(in_dir, out_dir, recurse.get(), want_vms.get(), want_csv.get()),
                    daemon=True).start()
                dlg.after(100, poll)

            run_btn = ttk.Button(frm, text="Run Bulk Correction", command=run,
                                 style="Accent.TButton")
            run_btn.grid(row=5, column=0, columnspan=3, pady=10)

    root = tk.Tk()
    style = ttk.Style()
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass
    style.configure("Accent.TButton", foreground="white", background="#005fb8",
                    font=("Helvetica", 9, "bold"))
    XPSViewerApp(root)
    root.mainloop()


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description="Transmission-correct CasaXPS VAMAS files (GUI or bulk CLI).")
    p.add_argument("--bulk", metavar="INPUT_DIR",
                   help="Correct every .vms file in this directory (headless).")
    p.add_argument("--out", metavar="OUTPUT_DIR",
                   help="Output directory (default: INPUT_DIR/corrected).")
    p.add_argument("--recurse", action="store_true",
                   help="Also process .vms files in subdirectories.")
    p.add_argument("--no-vms", action="store_true", help="Skip corrected VMS output.")
    p.add_argument("--no-csv", action="store_true", help="Skip comparison CSV output.")
    args = p.parse_args(argv)

    if args.bulk:
        in_dir = args.bulk
        if not os.path.isdir(in_dir):
            print(f"error: not a directory: {in_dir}", file=sys.stderr)
            return 2
        out_dir = args.out or os.path.join(in_dir, "corrected")
        summary = bulk_process_directory(
            in_dir, out_dir,
            recurse=args.recurse,
            write_vms=not args.no_vms,
            write_csv=not args.no_csv,
            log_cb=print)
        return 0 if summary["failed"] == 0 else 1

    try:
        launch_gui()
    except ImportError as exc:
        print(
            f"error: the GUI needs matplotlib and tkinter ({exc}).\n"
            "  - install GUI extras:  pip install 'harwellxps-vamas-correct[gui]'\n"
            "  - on minimal Linux you may also need Tk:  sudo apt install python3-tk\n"
            "  - or run headless:     harwell-xps --bulk INPUT_DIR",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
