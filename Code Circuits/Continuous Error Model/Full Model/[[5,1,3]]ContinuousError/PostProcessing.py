from qiskit import *
import numpy as np
import matplotlib.pyplot as plt

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, transpile
from qiskit.visualization import (
    plot_bloch_multivector,
    plot_histogram,
    latex,
    circuit_drawer,
)

from qiskit_aer import AerSimulator
from qiskit_aer.noise import (
    NoiseModel,
    QuantumError,
    ReadoutError,
    pauli_error,
    depolarizing_error,
    thermal_relaxation_error,
)

from qiskit.circuit import QuantumCircuit, Qubit, Clbit
import time

# ---------------------------------------------------------------------
# Logical codeword supports for the [[5,1,3]] code
# ---------------------------------------------------------------------

# Computational-basis states appearing in the logical |0_L> codeword.
L0 = [
    '00000','10010','01001','10100',
    '01010','11011','00110','11000',
    '11101','00011','11110','01111',
    '10001','01100','10111','00101'
]

# Computational-basis states appearing in the logical |1_L> codeword.
L1 = [
    '11111','01101','10110','01011',
    '10101','00100','11001','00111',
    '00010','11100','00001','10000',
    '01110','10011','01000','11010'
]


# ---------------------------------------------------------------------
# Logical decoding
# ---------------------------------------------------------------------

def post_processing_decoding(dic):
    """
    Decode measured physical states into logical states.

    Supports one or two encoded [[5,1,3]] logical qubits.

    Parameters
    ----------
    dic : dict
        Raw Qiskit count dictionary.

    Returns
    -------
    dict
        Logical count dictionary.
    """

    # --------------------------------------------------------------
    # One encoded logical qubit
    # --------------------------------------------------------------
    if len(next(iter(dic))) < 10:

        dic_out = {'0': 0, '1': 0}

        for item, count in dic.items():

            # Reverse the Qiskit bit string to recover the ordering
            # of the physical data qubits.
            iteminv = item[::-1]

            logical_bits = iteminv[:5][::-1]

            if logical_bits in L0:
                dic_out['0'] += count

            elif logical_bits in L1:
                dic_out['1'] += count

    # --------------------------------------------------------------
    # Two encoded logical qubits
    # --------------------------------------------------------------
    elif len(next(iter(dic))) == 10 or len(next(iter(dic))) == 14:

        dic_out = {
            '00': 0,
            '01': 0,
            '10': 0,
            '11': 0
        }

        for key, value in dic.items():

            keyinv = key[::-1]

            # Extract the two encoded blocks.
            first_5_key = keyinv[:5][::-1]
            last_5_key = keyinv[5:10][::-1]

            # Decode each logical block independently.
            if first_5_key in L0 and last_5_key in L0:
                dic_out['00'] += value

            elif first_5_key in L0 and last_5_key in L1:
                dic_out['01'] += value

            elif first_5_key in L1 and last_5_key in L0:
                dic_out['10'] += value

            elif first_5_key in L1 and last_5_key in L1:
                dic_out['11'] += value

    return dic_out


# ---------------------------------------------------------------------
# Ideal circuit simulation
# ---------------------------------------------------------------------

def QuantumSimulator(qc, shots=100):
    """
    Execute an ideal quantum circuit.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to simulate.

    shots : int
        Number of measurement shots.

    Returns
    -------
    dict
        Measurement counts.
    """

    simulator = AerSimulator()

    job = simulator.run(
        qc,
        shots=shots
    )

    result = job.result()

    measurement_counts = result.get_counts()

    return measurement_counts


# ---------------------------------------------------------------------
# Dictionary utilities
# ---------------------------------------------------------------------

def Sum_dict(dict1, dict2):
    """
    Merge two count dictionaries by summing common entries.

    Parameters
    ----------
    dict1, dict2 : dict

    Returns
    -------
    dict
        Combined count dictionary.
    """

    merged_dict = dict1.copy()

    for key, value in dict2.items():

        if key in merged_dict:
            merged_dict[key] += value

        else:
            merged_dict[key] = value

    return merged_dict


# ---------------------------------------------------------------------
# Repeated ideal simulations
# ---------------------------------------------------------------------

def runcircuitnoisy(qc, shots):
    """
    Repeat the circuit several times and accumulate the counts.

    Although the name contains 'noisy', this function currently performs
    repeated ideal simulations. Each iteration executes 100 shots and the
    resulting count dictionaries are accumulated.

    Parameters
    ----------
    qc : QuantumCircuit

    shots : int
        Number of repetitions.

    Returns
    -------
    dict
        Total accumulated counts.
    """

    counts = {}

    for _ in range(shots):

        # Run one batch of 100 ideal shots.
        r = QuantumSimulator(qc, 100)

        # Add the new counts to the accumulated dictionary.
        counts = Sum_dict(counts, r)

    return counts


# ---------------------------------------------------------------------
# Repeat simulations and decode logical outcomes
# ---------------------------------------------------------------------

def runcircxtimes(qc, shots, times):
    """
    Repeat the simulation multiple times and decode each result
    into logical measurement outcomes.

    Parameters
    ----------
    qc : QuantumCircuit

    shots : int
        Number of repetitions used inside runcircuitnoisy().

    times : int
        Number of independent simulations.

    Returns
    -------
    list
        List of decoded logical count dictionaries.
    """

    prob = []

    for _ in range(int(times)):

        # Run one complete simulation.
        counts = runcircuitnoisy(qc, int(shots))

        # Decode physical outcomes into logical outcomes.
        countsL = post_processing_decoding(counts)

        prob.append(countsL)

    return prob
