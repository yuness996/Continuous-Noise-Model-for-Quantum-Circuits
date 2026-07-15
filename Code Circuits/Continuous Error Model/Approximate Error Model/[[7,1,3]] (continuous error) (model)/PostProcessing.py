from qiskit import *
import numpy as np

from qiskit import (
    QuantumRegister,
    ClassicalRegister,
    QuantumCircuit,
    transpile,
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

# ---------------------------------------------------------------------
# Logical codeword supports for the [[7,1,3]] Steane code
# ---------------------------------------------------------------------

# Computational-basis states appearing in the logical |0_L> codeword.
L0 = [
    '0000000','1010101','0110011','1100110',
    '0001111','1011010','0111100','1101001'
]

# Computational-basis states appearing in the logical |1_L> codeword.
L1 = [
    '1111111','0101010','1001100','0011001',
    '1110000','0100101','1000011','0010110'
]


# ---------------------------------------------------------------------
# Logical decoding
# ---------------------------------------------------------------------

def post_processing_decoding(dic):
    """
    Decode physical measurement outcomes into logical outcomes.

    Supports one or two encoded [[7,1,3]] logical qubits.

    Parameters
    ----------
    dic : dict
        Raw Qiskit measurement counts.

    Returns
    -------
    dict
        Dictionary containing logical measurement counts.
    """

    # --------------------------------------------------------------
    # One encoded logical qubit
    # --------------------------------------------------------------
    if len(next(iter(dic))) < 14:

        dic_out = {'0': 0, '1': 0}

        for item, count in dic.items():

            # Recover the seven measured data bits.
            itemL = item[::-1][:7][::-1]

            # Decode the measured codeword.
            if itemL in L0:
                dic_out['0'] += count

            else:
                dic_out['1'] += count

    # --------------------------------------------------------------
    # Two encoded logical qubits
    # --------------------------------------------------------------
    elif len(next(iter(dic))) == 14 or len(next(iter(dic))) == 20:

        dic_out = {
            '00': 0,
            '01': 0,
            '10': 0,
            '11': 0
        }

        for key, value in dic.items():

            # The bit ordering depends on whether syndrome bits
            # are present in the measurement string.
            if len(next(iter(dic))) == 20:

                state_qubit_1 = key[::-1][:7][::-1]
                state_qubit_2 = key[::-1][13:20][::-1]

            elif len(next(iter(dic))) == 14:

                state_qubit_1 = key[7:]
                state_qubit_2 = key[:7]

            # Decode the two logical qubits.
            if state_qubit_1 in L0 and state_qubit_2 in L0:

                dic_out['00'] += value

            elif state_qubit_1 in L1 and state_qubit_2 in L1:

                dic_out['11'] += value

            elif (state_qubit_1 in L0) or (state_qubit_2 in L1):

                dic_out['01'] += value

            elif (state_qubit_1 in L1) or (state_qubit_2 in L0):

                dic_out['10'] += value

            # If one or both measured blocks do not belong to the
            # logical supports, split the counts equally between the
            # two mixed logical states.
            else:

                dic_out['01'] += int(value / 2)
                dic_out['10'] += value - int(value / 2)

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
        Raw measurement counts.
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
    Merge two count dictionaries by summing matching entries.

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
# Repeated ideal simulation
# ---------------------------------------------------------------------

def runcircuitnoisy(qc, shots):
    """
    Execute the circuit repeatedly and accumulate measurement counts.

    Despite its name, this function performs repeated ideal simulations.
    Each repetition consists of 100 measurement shots.

    Parameters
    ----------
    qc : QuantumCircuit

    shots : int
        Number of repetitions.

    Returns
    -------
    dict
        Accumulated measurement counts.
    """

    counts = {}

    for _ in range(shots):

        # Execute one batch of 100 shots.
        result = QuantumSimulator(qc, 100)

        # Accumulate the measurement counts.
        counts = Sum_dict(counts, result)

    return counts


# ---------------------------------------------------------------------
# Repeat simulations and decode logical outcomes
# ---------------------------------------------------------------------

def runcircxtimes(qc, shots, times):
    """
    Repeat the circuit simulation and decode the results into logical
    measurement outcomes.

    Parameters
    ----------
    qc : QuantumCircuit

    shots : int
        Number of repetitions used by runcircuitnoisy().

    times : int
        Number of independent simulations.

    Returns
    -------
    list
        List of decoded logical count dictionaries.
    """

    prob = []

    for _ in range(int(times)):

        # Execute one complete simulation.
        counts = runcircuitnoisy(qc, int(shots))

        # Decode physical codewords into logical states.
        countsL = post_processing_decoding(counts)

        # Uncomment if normalized probabilities are required.
        # countsL = normalize(countsL, int(shots) * 100)

        prob.append(countsL)

    return prob
