"""Run the noisy [[5,1,3]] encoded-qubit circuit over a range of noise strengths."""

import os

import numpy as np

from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

# Circuit blocks for encoding, logical gates, syndrome extraction, and recovery.
from Oracles import *

# Simulation and logical post-processing functions.
from PostProcessing import *

# ---------------------------------------------------------------------
# Single encoded-qubit simulation
# ---------------------------------------------------------------------

def single_encoded_qubit(k, shots, times, n_h, correction):
    """
    Build and simulate a noisy [[5,1,3]] encoded-qubit circuit.

    The circuit encodes one logical qubit, applies ten noisy logical
    Hadamard circuits, and measures the five data qubits.

    Error correction can be enabled by allocating four syndrome ancillas.
    However, the syndrome-extraction and recovery blocks are currently
    commented out, so setting ``correction=True`` only increases the
    register size.

    Parameters
    ----------
    k : float
        Strength of the continuous coherent noise used by ``encode()``
        and ``Logical_Hadamard()``.

    shots : int
        Number of measurement shots used in each simulation batch.

    times : int
        Number of independent noisy circuit instances.

    correction : bool
        Whether four syndrome ancillas should be allocated.

    Returns
    -------
    object
        Result returned by ``runcircxtimes()``.
    """

    # The [[5,1,3]] code uses five data qubits.
    # Four additional ancillas are required for syndrome extraction.
    if correction:
        nbq = 4
    else:
        nbq = 0

    # Create quantum and classical registers.
    qq = QuantumRegister(5 + nbq)
    cc = ClassicalRegister(5 + nbq)

    qc = QuantumCircuit(qq, cc)

    # Add the noisy [[5,1,3]] encoding circuit.
    qc.compose(
        encode(k),
        range(5),
        inplace=True,
    )

    # Apply ten noisy logical Hadamard circuits.
    for _ in range(n_h):
        qc.compose(Logical_Hadamard(k), range(5), inplace=True)

    # The following block would perform one round of syndrome extraction and correction after the logical Hadamards.
        
    if correction:
        qc.compose(detection(k), range(9), inplace=True)
        qc = decode_correction(qc, qq, [cc[5], cc[6], cc[7], cc[8]])

    # Measure the five physical data qubits.
    for i in range(5):
        qc.measure(qq[i], cc[i])

    # Repeat the simulation and compare the decoded output with the
    # expected logical distribution.
    return runcircxtimes(qc, shots, times)


# ---------------------------------------------------------------------
# Sweep over noise strengths
# ---------------------------------------------------------------------

def Run1Qbt_AllNoise(a, b, num, shots, times, n_h, correction):
    """
    Simulate the encoded-qubit circuit over a logarithmic noise grid.

    Parameters
    ----------
    a, b : float
        Lower and upper exponents of the noise range. The tested values
        run from 10**a to 10**b.

    num : int
        Number of noise strengths.

    shots : int
        Number of shots per simulation batch.

    times : int
        Number of independent noisy circuit instances per noise value.

    correction : bool
        Whether syndrome ancillas should be allocated.

    Returns
    -------
    list
        Simulation results for all tested noise strengths.
    """

    # Generate logarithmically spaced continuous-noise strengths.
    noise_strengths = np.logspace(
        a,
        b,
        num=num,
    )

    results = []

    for k in noise_strengths:
        results.append(single_encoded_qubit(k, shots=shots, times=times, n_h=n_h, correction=correction))

    return results


# ---------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------

def save_table(table, filename):
    """
    Append simulation results to a text file.

    Each row is written on a separate line, with entries separated
    by spaces.

    Parameters
    ----------
    table : iterable
        Simulation results to save.

    filename : str
        Output file name.
    """

    with open(filename, "a", encoding="utf-8") as file:

        # Optional header:
        # file.write("Shots = 1000\n")

        for row in table:
            file.write(" ".join(map(str, row)) + "\n")


# ---------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------

if __name__ == "__main__":
    results = Run1Qbt_AllNoise(a=-6, b=-1, num=20, shots=1000, times=10, n_h=10, correction=True)

    # Append the results to the output file.
    save_table(results, "test.txt")

