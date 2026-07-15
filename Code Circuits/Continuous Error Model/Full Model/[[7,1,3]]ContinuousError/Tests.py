"""Build and scan the noisy [[7,1,3]] encoded-qubit circuit."""

import numpy as np
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister

from Oracles import decode_correction, detection, encode
from PostProcessing import *


def single_encoded_qubit(k, shots, times, correction):
    """
    Build one noisy [[7,1,3]] encoded-qubit circuit.

    The circuit encodes one logical qubit, optionally performs one round
    of syndrome extraction and correction, and measures the seven data
    qubits.

    ``shots`` and ``times`` are retained for compatibility but are not
    used because this function currently returns the circuit itself.
    """
    number_of_ancillas = 6 if correction else 0

    qreg = QuantumRegister(7 + number_of_ancillas, "q")
    creg = ClassicalRegister(7 + number_of_ancillas, "c")
    circuit = QuantumCircuit(qreg, creg)

    # Add the noisy encoding circuit.
    circuit.compose(encode(k), range(7), inplace=True)

    # Apply one noisy syndrome-extraction and recovery cycle.
    if correction:
        circuit.compose(detection(k), range(13), inplace=True)
        circuit = decode_correction(circuit, qreg, creg)

    # Measure the seven physical data qubits.
    for qubit in range(7):
        circuit.measure(qreg[qubit], creg[qubit])

    # Repeat the simulation and compare the decoded output with the
    # expected logical distribution.
    return runcircxtimes(circuit, shots, times)


def Run1Qbt_AllNoise(a, b, num, shots, times, correction):
    """Build one circuit for each noise strength in a logarithmic grid."""
    noise_strengths = np.logspace(a, b, num=num)

    return [
        single_encoded_qubit(
            k,
            shots=shots,
            times=times,
            correction=correction,
        )
        for k in noise_strengths
    ]


def save_table(table, filename):
    """
    Append rows of numerical results to a text file.

    This function expects each element of ``table`` to be iterable.
    It cannot directly save the QuantumCircuit objects currently returned
    by ``Run1Qbt_AllNoise``.
    """
    with open(filename, "a", encoding="utf-8") as file:
        file.write("Shots = 1000\n")

        for row in table:
            file.write(" ".join(map(str, row)) + "\n")


if __name__ == "__main__":
    results = Run1Qbt_AllNoise(
        a=-6,
        b=-1,
        num=20,
        shots=10,
        times=2,
        correction=True,
    )
    # Append the results to the output file.
    save_table(results, "test.txt")