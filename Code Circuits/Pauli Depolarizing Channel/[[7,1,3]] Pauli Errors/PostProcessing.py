import numpy as np

from itertools import product

from qiskit import transpile
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, pauli_error


# ---------------------------------------------------------------------
# Logical codeword supports for the [[7,1,3]] code
# ---------------------------------------------------------------------

L0 = [
    "0000000", "1010101", "0110011", "1100110",
    "0001111", "1011010", "0111100", "1101001",
]

L1 = [
    "1111111", "0101010", "1001100", "0011001",
    "1110000", "0100101", "1000011", "0010110",
]


# ---------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------

def post_processing_decoding(dic):
    """
    Decode measured physical states into logical states.

    This version is for one encoded [[7,1,3]] qubit.

    Assumption
    ----------
    The seven data-qubit measurement bits are the rightmost seven bits
    in the Qiskit count key.
    """
    if len(next(iter(dic))) < 14:
        dic_out = {"0": 0, "1": 0}

        for item, count in dic.items():
            iteminv = item[::-1]
            itemL = iteminv[:7][::-1]

            if itemL in L0:
                dic_out["0"] += count

            else:
                dic_out["1"] += count

    return dic_out


# ---------------------------------------------------------------------
# Noisy simulation
# ---------------------------------------------------------------------

def runcircuitnoisy(qc, err, p_err, shots):
    """
    Run a quantum circuit qc in a Pauli noise channel.

    err : list[str]
        Type of Pauli errors, for example ["X", "Y", "Z"].

    p_err : list[float]
        Probability of each Pauli error.

    shots : int
        Number of shots.
    """
    arr = []

    for i in range(len(err)):
        arr.append((err[i], p_err[i]))

    arr.append(("I", 1 - np.sum(p_err)))

    error = pauli_error(arr)

    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(error, ["u1", "u2", "u3"])

    sim_noise = AerSimulator(noise_model=noise_model)

    qc_noise = transpile(
        qc,
        sim_noise,
        optimization_level=0,
    )

    result = sim_noise.run(qc_noise, shots=int(shots)).result()
    counts = result.get_counts(0)

    return counts


# ---------------------------------------------------------------------
# Repeated noisy runs
# ---------------------------------------------------------------------

def runcircxtimes(qc, err, p_err, shots, times):
    """
    Run qc several times and return decoded logical count dictionaries.
    """
    prob = []

    for i in range(int(times)):
        counts = runcircuitnoisy(qc, err, p_err, int(shots))
        countsL = post_processing_decoding(counts)
        prob.append(countsL)

    return prob
