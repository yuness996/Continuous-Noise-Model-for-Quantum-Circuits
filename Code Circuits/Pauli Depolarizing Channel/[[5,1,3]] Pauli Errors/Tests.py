import numpy as np

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit

# Contains the circuit blocks used below, such as:
# encode(), Logical_Hadamard(), and detection().
from Oracles import *

# Contains the simulation and correction tools, such as:
# decode_correction() and runcircxtimes().
from PostProcessing import *


def single_encoded_qubit(p_err, shots, times, correction, n_h):
    """
    Build and simulate a [[5,1,3]] encoded-qubit circuit under Pauli noise.

    The circuit:
    1. Encodes one logical qubit into five physical qubits.
    2. Applies the logical Hadamard circuit n_h times.
    3. Optionally performs syndrome extraction and error correction.
    4. Measures the five data qubits.
    5. Repeats the noisy simulation several times.

    Parameters
    ----------
    p_err : float
        Probability assigned to each Pauli error X, Y, and Z.
        The total depolarizing error probability is therefore 3 * p_err.

    shots : int
        Number of measurement shots used for each noisy circuit instance.

    times : int
        Number of independent noisy circuit instances to simulate.

    correction : bool
        If True, add four syndrome qubits and perform error correction.

    n_h : int
        Number of logical Hadamard circuits applied after encoding.

    Returns
    -------
    Simulation result returned by runcircxtimes().
    """

    # The [[5,1,3]] code uses five data qubits.
    # Four additional ancilla qubits are required for syndrome extraction.
    if correction:
        fn = 4
    else:
        fn = 0

    # Create quantum and classical registers.
    # Their size is 5 without correction and 9 with correction.
    qq = QuantumRegister(5 + fn, "q")
    cc = ClassicalRegister(5 + fn, "c")
    qc = QuantumCircuit(qq, cc)

    # Encode the initial logical qubit into the five-qubit code space.
    qc.append(encode(), range(5))

    # Apply the logical Hadamard circuit repeatedly.
    # Barriers separate consecutive logical operations in the circuit diagram.
    for _ in range(n_h):
        qc.append(Logical_Hadamard(), range(5))
        qc.barrier()

    if correction:
        # Perform syndrome extraction using the five data qubits
        # and four syndrome ancillas.
        qc.append(detection(), range(9))

        # Measure the syndrome qubits and apply the correction associated
        # with the measured syndrome.
        qc = decode_correction(
            qc,
            qq,
            [cc[5], cc[6], cc[7], cc[8]]
        )

    # Measure the five physical data qubits.
    for i in range(5):
        qc.measure(qq[i], cc[i])

    # Run the circuit under a symmetric Pauli channel.
    # X, Y, and Z errors are each applied with probability p_err.
    return runcircxtimes(
        qc,
        ["X", "Y", "Z"],
        [p_err] * 3,
        shots,
        times
    )


def Run_AllNoise(a, b, num, shots, times, correction, n_h):
    """
    Simulate the encoded-qubit circuit over a logarithmic noise grid.

    Parameters
    ----------
    a, b : float
        Lower and upper exponents of the total Pauli error-rate range.
        The tested values are generated from 10**a to 10**b.

    num : int
        Number of noise values in the logarithmic grid.

    shots : int
        Number of measurement shots per circuit instance.

    times : int
        Number of independent noisy circuit instances per noise value.

    correction : bool
        Whether syndrome extraction and correction are applied.

    n_h : int
        Number of logical Hadamard circuits.

    Returns
    -------
    list
        Simulation results for all tested physical error rates.
    """

    # Generate total depolarizing error probabilities on a logarithmic grid.
    p_errors = np.logspace(a, b, num=num)

    results = []

    for p_total in p_errors:
        # Divide the total Pauli error probability equally among X, Y, and Z.
        results.append(
            single_encoded_qubit(
                p_err=p_total / 3,
                shots=shots,
                times=times,
                correction=correction,
                n_h=n_h
            )
        )

    return results


def save_table(table, filename):
    """
    Append a table of results to a text file.

    Each row is written on a separate line, with values separated by spaces.
    """

    with open(filename, "a", encoding="utf-8") as file:
        for row in table:
            file.write(" ".join(map(str, row)) + "\n")


if __name__ == "__main__":
    # Run the circuit for 20 total Pauli error rates between 10^-6 and 10^-1.
    #
    # In this example:
    # - no error correction is applied,
    # - no logical Hadamard is applied,
    # - each circuit uses 1000 shots,
    # - each noise value is simulated using 10 independent instances.
    results = Run_AllNoise(
        a=-6,
        b=-1,
        num=20,
        shots=1000,
        times=10,
        correction=False,
        n_h=0
    )

    # Append the simulation results to the output file.
    save_table(results, "test.txt")
