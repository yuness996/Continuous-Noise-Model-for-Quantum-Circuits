"""Read repeated GSA simulation results and compute their column-wise mean.

The input file may contain metadata lines beginning with ``Shots`` followed by
one or more rows of whitespace-separated floating-point values. Each numerical
row is expected to contain the same number of points from a noise sweep.
"""

from pathlib import Path

import numpy as np


# Number of noise-strength points written by the default GSA simulations.
DEFAULT_NUM_POINTS = 20


def read_mean_array(filename: str | Path, expected_length: int = DEFAULT_NUM_POINTS) -> np.ndarray:
    """Load numerical result rows and return their column-wise average.

    Parameters
    ----------
    filename:
        Path to the text file containing simulation outputs.
    expected_length:
        Required number of floating-point values in each numerical row.

    Returns
    -------
    numpy.ndarray
        One-dimensional array containing the mean of each column.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist.
    ValueError
        If ``expected_length`` is not positive, a numerical row has the wrong
        number of entries, a value cannot be converted to ``float``, or the
        file contains no numerical data rows.
    """
    if expected_length <= 0:
        raise ValueError("expected_length must be a positive integer.")

    rows: list[list[float]] = []
    input_path = Path(filename)

    with input_path.open("r", encoding="utf-8") as file:
        for line_number, line in enumerate(file, start=1):
            stripped_line = line.strip()

            # Ignore blank lines and metadata written by the simulation script.
            if not stripped_line or stripped_line.startswith("Shots"):
                continue

            try:
                values = [float(value) for value in stripped_line.split()]
            except ValueError as exc:
                raise ValueError(
                    f"Line {line_number} in '{input_path}' contains a non-numeric value."
                ) from exc

            # Every row must describe the same noise sweep.
            if len(values) != expected_length:
                raise ValueError(
                    f"Line {line_number} in '{input_path}' contains {len(values)} "
                    f"values; expected {expected_length}."
                )

            rows.append(values)

    if not rows:
        raise ValueError(f"No numerical data rows were found in '{input_path}'.")

    # Axis 0 averages repeated runs point by point across the noise sweep.
    return np.mean(np.asarray(rows, dtype=float), axis=0)


def main() -> None:
    """Read the default seven-qubit result file and print the mean as a list."""
    mean_values = read_mean_array("GSA7qubits_cont.txt")
    print(mean_values.tolist())


if __name__ == "__main__":
    main()