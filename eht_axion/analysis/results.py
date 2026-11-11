"""Result files of the limit stages: numeric arrays in results.npz, everything else in report.json.

A stage writes one report (identity, summary table, notes) plus named sections of records,
for example ``observations`` {"1d": {...}, "2d": {...}} or ``cases`` {"observed-cor-1d": {...}}.
Every nonempty rectangular numeric list or array of a record is stored once, in results.npz under
"<section>/<name>/<field>"; the remaining fields stay in report.json under report[section][name].
read_results returns the same records with the arrays as NumPy arrays.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPORT, ARRAYS = "report.json", "results.npz"


def _array(value):
    if isinstance(value, np.ndarray):
        return value
    if not isinstance(value, (list, tuple)) or not value:
        return None
    try:
        array = np.asarray(value)
    except ValueError:  # ragged
        return None
    return array if array.dtype.kind in "biuf" and array.size else None


def write_results(directory, report, sections):
    """Write report.json and results.npz of one stage directory."""
    directory = Path(directory)
    text, arrays = dict(report), {}
    for section, records in sections.items():
        if section in text:
            raise ValueError(f"report field and section collide: {section}")
        text[section] = {}
        for name, record in records.items():
            kept = {}
            for field, value in record.items():
                array = _array(value)
                if array is None:
                    kept[field] = value
                else:
                    arrays[f"{section}/{name}/{field}"] = array
            text[section][name] = kept
    np.savez(directory / ARRAYS, **arrays)
    (directory / REPORT).write_text(json.dumps(text, indent=2, allow_nan=False) + "\n")
    return [directory / REPORT, directory / ARRAYS]


def read_results(directory, sections):
    """Return (report without the sections, {section: {name: record}}) of one stage directory."""
    directory = Path(directory)
    report = json.loads((directory / REPORT).read_text())
    records = {section: {name: dict(fields) for name, fields in report.pop(section).items()}
               for section in sections}
    with np.load(directory / ARRAYS, allow_pickle=False) as stored:
        for key in stored.files:
            section, name, field = key.split("/")
            if section in records:
                records[section][name][field] = stored[key]
    return report, records
