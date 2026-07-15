"""
Utilities for reading one-encoded-qubit simulation result files.

Expected file format:

    Shots = 1000
    {'0': 100000, '1': 0} {'0': 99999, '1': 1} ...
    {'0': 99950, '1': 50} {'0': 99960, '1': 40} ...
    ...

Each non-header line corresponds to one noise value.
Each line contains several count dictionaries, usually one per repeated run.

Important:
    The probability is normalized using sum(counts.values()),
    not the "Shots =" header.
"""

import ast
import re
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------

def parse_dicts_from_line(line):
    """
    Extract all count dictionaries from one line.

    Example
    -------
    "{'0': 10, '1': 2} {'0': 9, '1': 3}"

    becomes

    [{'0': 10, '1': 2}, {'0': 9, '1': 3}]
    """
    dict_strings = re.findall(r"\{[^{}]*\}", line)

    dictionaries = []

    for dict_string in dict_strings:
        dictionaries.append(ast.literal_eval(dict_string))

    return dictionaries


def read_result_rows(filename, header_prefix="Shots"):
    """
    Read result rows from a file.

    Lines starting with "Shots" are skipped.
    Empty lines are skipped.

    Returns
    -------
    rows : list[list[dict]]
        rows[noise_index][repeat_index] = count dictionary

    n_headers : int
        Number of detected "Shots =" blocks.
    """
    filename = Path(filename)

    rows = []
    n_headers = 0

    with filename.open("r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(header_prefix):
                n_headers += 1
                continue

            parsed = parse_dicts_from_line(line)

            if parsed:
                rows.append(parsed)

    return rows, n_headers


# ---------------------------------------------------------------------
# Grouping and merging parallel jobs
# ---------------------------------------------------------------------

def group_rows_by_sweep(rows, sweep_size):
    """
    Group rows into independent jobs.

    If one job contains `sweep_size` noise values, then the rows are grouped as

        rows[0:sweep_size],
        rows[sweep_size:2*sweep_size],
        ...

    Returns
    -------
    grouped : list
        grouped[job][noise_index][repeat_index]
    """
    if len(rows) % sweep_size != 0:
        raise ValueError(
            f"Number of data rows ({len(rows)}) is not divisible by "
            f"sweep_size ({sweep_size})."
        )

    grouped = []

    for start in range(0, len(rows), sweep_size):
        grouped.append(rows[start:start + sweep_size])

    return grouped


def add_count_dictionaries(dict_1, dict_2):
    """
    Add two count dictionaries.

    Missing keys are treated as zero.
    """
    keys = set(dict_1) | set(dict_2)

    return {
        key: dict_1.get(key, 0) + dict_2.get(key, 0)
        for key in keys
    }


def merge_parallel_jobs(grouped_rows):
    """
    Merge several jobs by adding count dictionaries.

    Input shape:
        grouped_rows[job][noise_index][repeat_index]

    Output shape:
        merged[noise_index][repeat_index]
    """
    if not grouped_rows:
        return []

    n_noise = len(grouped_rows[0])
    n_repeats = len(grouped_rows[0][0])

    merged = [
        [dict() for _ in range(n_repeats)]
        for _ in range(n_noise)
    ]

    for job in grouped_rows:
        if len(job) != n_noise:
            raise ValueError("All jobs must have the same number of noise rows.")

        for noise_index in range(n_noise):
            if len(job[noise_index]) != n_repeats:
                raise ValueError(
                    "All noise rows must have the same number of repeated dictionaries."
                )

            for repeat_index in range(n_repeats):
                merged[noise_index][repeat_index] = add_count_dictionaries(
                    merged[noise_index][repeat_index],
                    job[noise_index][repeat_index],
                )

    return merged


# ---------------------------------------------------------------------
# One-encoded-qubit observable
# ---------------------------------------------------------------------

def logical_error_1q(counts):
    """
    Logical error probability for one encoded qubit.

    For one encoded qubit, the logical error is

        P_L = counts['1'] / total_counts

    where

        total_counts = counts['0'] + counts['1'].

    This avoids probabilities larger than one when the file contains
    dictionaries such as {'0': 100000, '1': 0}.
    """
    total_counts = sum(counts.values())

    if total_counts <= 0:
        raise ValueError(f"Invalid count dictionary with zero total: {counts}")

    n_error = counts.get("1", 0) + counts.get(1, 0)

    return n_error / total_counts


def mean_and_error_bar(values):
    """
    Return the sample mean and the standard error of the mean.
    """
    values = np.asarray(values, dtype=float)

    mean = float(np.mean(values))

    if len(values) <= 1:
        return mean, 0.0

    std = float(np.std(values, ddof=1))
    error_bar = std / np.sqrt(len(values))

    return mean, error_bar


# ---------------------------------------------------------------------
# Main post-processing
# ---------------------------------------------------------------------

def compute_error_curve(
    merged_rows,
    times=None,
    prepend_zero=True,
):
    """
    Compute the logical error curve for one encoded qubit.

    Parameters
    ----------
    merged_rows : list
        Output of merge_parallel_jobs(...), with shape

            merged_rows[noise_index][repeat_index]

    times : int or None
        Number of repeated dictionaries to use per noise value.
        If None, all dictionaries in each row are used.

    prepend_zero : bool
        If True, prepend zero to the probability and error-bar arrays.

    Returns
    -------
    probabilities : np.ndarray
        Mean logical error probability for each noise value.

    error_bars : np.ndarray
        Standard error bar for each noise value.
    """
    probabilities = []
    error_bars = []

    for noise_row in merged_rows:
        if times is None:
            selected_counts = noise_row
        else:
            if times > len(noise_row):
                raise ValueError(
                    f"times={times} but this row only contains "
                    f"{len(noise_row)} dictionaries."
                )

            selected_counts = noise_row[:times]

        values = [
            logical_error_1q(counts)
            for counts in selected_counts
        ]

        mean, err = mean_and_error_bar(values)

        probabilities.append(mean)
        error_bars.append(err)

    probabilities = np.asarray(probabilities, dtype=float)
    error_bars = np.asarray(error_bars, dtype=float)

    if prepend_zero:
        probabilities = np.concatenate(([0.0], probabilities))
        error_bars = np.concatenate(([0.0], error_bars))

    return probabilities, error_bars


def load_error_curve(
    filename,
    sweep_size,
    times=None,
    prepend_zero=True,
):
    """
    Load a one-encoded-qubit result file and compute the logical error curve.

    Parameters
    ----------
    filename : str
        Path to the result file.

    sweep_size : int
        Number of noise values in one sweep.
        For your files this is usually 20.

    times : int or None
        Number of repeated dictionaries per noise value.
        For your files this is usually 10.
        If None, all dictionaries in each row are used.

    prepend_zero : bool
        If True, prepend zero to the output arrays.

    Returns
    -------
    probabilities : np.ndarray
        Mean logical error probabilities.

    error_bars : np.ndarray
        Standard error bars.
    """
    rows, n_headers = read_result_rows(filename)

    grouped = group_rows_by_sweep(
        rows=rows,
        sweep_size=sweep_size,
    )

    merged = merge_parallel_jobs(grouped)

    probabilities, error_bars = compute_error_curve(
        merged_rows=merged,
        times=times,
        prepend_zero=prepend_zero,
    )

    return probabilities, error_bars


# ---------------------------------------------------------------------
# Example
# ---------------------------------------------------------------------

fn = "test.txt"

curve, error_bar = load_error_curve(
    filename=fn,
    sweep_size=20,
    times=10,
    prepend_zero=True,
)

print("curve = [" + ",".join(map(str, curve)) + "]")
print("error_bar = [" + ",".join(map(str, error_bar)) + "]")