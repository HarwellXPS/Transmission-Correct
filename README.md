name: CI

on:
  push:
    branches: [main, master]
  pull_request:

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        python-version: ["3.8", "3.9", "3.10", "3.11", "3.12"]

    steps:
      - uses: actions/checkout@v4

      - name: Set up Python ${{ matrix.python-version }}
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install package (core + test extras)
        run: |
          python -m pip install --upgrade pip
          pip install -e ".[test]"

      - name: Byte-compile
        run: python -m py_compile xps_transmission_tool.py

      - name: Run tests
        run: pytest -q
