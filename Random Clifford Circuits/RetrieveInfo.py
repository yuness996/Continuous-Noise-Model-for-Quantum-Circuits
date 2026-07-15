import numpy as np


def load_piled_matrices(filename, matrix_shape=(10, 10), ddof=0):
    matrices = []
    current_matrix = []

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            row = [float(x) for x in line.split()]
            current_matrix.append(row)

            if len(current_matrix) == matrix_shape[0]:
                matrix = np.array(current_matrix, dtype=float)

                if matrix.shape != matrix_shape:
                    raise ValueError(
                        f"Expected matrix shape {matrix_shape}, got {matrix.shape}"
                    )

                matrices.append(matrix)
                current_matrix = []

    if current_matrix:
        raise ValueError(
            f"Incomplete matrix: found {len(current_matrix)} rows "
            f"instead of {matrix_shape[0]}."
        )

    matrices = np.array(matrices)

    if matrices.size == 0:
        raise ValueError("No matrices were found in the file.")

    mean_matrix = np.mean(matrices, axis=0)
    variance_matrix = np.var(matrices, axis=0, ddof=ddof)

    return matrices, mean_matrix, variance_matrix


filename = "RandomCircuits10qbts_new (copy).txt"

matrices, mean_matrix, variance_matrix = load_piled_matrices(
    filename,
    matrix_shape=(10, 10),
    ddof=0
)

print("Number of matrices loaded:", len(matrices))

print("Mean matrix:")
print(np.array2string(
    mean_matrix,
    separator=", ",
    formatter={"float_kind": lambda x: f"{x:.6e}"}
))

print("Variance matrix:")
print(np.array2string(
    variance_matrix,
    separator=", ",
    formatter={"float_kind": lambda x: f"{x:.6e}"}
))