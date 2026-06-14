# Contributing

Thanks for your interest in improving this tool. Issues and pull requests are welcome.

## Reporting parsing problems

The VAMAS block parser uses pragmatic field-offset heuristics tuned to CasaXPS-style
exports. If a file returns **zero blocks**, the wrong point counts, or "block boundaries
could not be verified", the field layout likely differs from what the parser expects.

The fastest way to get support added is to open an issue with a **small, anonymised**
sample `.vms` file (a single region with a handful of points is enough) and a note on the
instrument/software that produced it.

## Development setup

```bash
git clone https://github.com/your-username/harwellxps-vamas-correct.git
cd harwellxps-vamas-correct
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[test,gui]"
```

## Running the tests

```bash
pytest
```

The test suite exercises the GUI-free core (parsing, correction, subset writing, bulk
processing) and needs only NumPy + pytest — no display required.

## Coding notes

- Keep the parsing/correction **core** free of GUI imports; `tkinter`/`matplotlib` are
  imported lazily inside `launch_gui()` so the CLI stays headless.
- Instrument-specific behaviour lives in constants near the top of
  `xps_transmission_tool.py` (`N_CORR_VARS`, `REGION_NAMES`, `SOURCE_LABELS`).
- Add or update a test in `tests/` for any behaviour change.
