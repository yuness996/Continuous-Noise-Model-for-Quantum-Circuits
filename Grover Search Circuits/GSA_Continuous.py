import random

import matplotlib.pyplot as plt
import numpy as np

from qiskit import (
    ClassicalRegister,
    QuantumCircuit,
    QuantumRegister,
    transpile,
)
from qiskit.circuit.library import UGate, ZGate
from qiskit.visualization import (
    circuit_drawer,
    latex,
    plot_bloch_multivector,
    plot_histogram,
)
from qiskit_aer import AerSimulator
from qiskit_aer.noise import (
    NoiseModel,
    QuantumError,
    ReadoutError,
    depolarizing_error,
    pauli_error,
    thermal_relaxation_error,
)


# ---------------------------------------------------------------------
# Continuous coherent-error model
# ---------------------------------------------------------------------

def r_error_gate(qc, theta, phi, qubit):
    """
    Apply the continuous coherent-error unitary to one qubit.

    The error is decomposed as

        Z Rx(2 theta) Rz(-2 phi) Z,

    where theta and phi are random rotation angles.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to which the error gate is appended.
    theta : float
        Rotation angle associated with the Rx component.
    phi : float
        Rotation angle associated with the Rz component.
    qubit : int
        Index of the affected qubit.

    Returns
    -------
    QuantumCircuit
        The modified circuit.
    """
    qc.z(qubit)
    qc.rx(2 * theta, qubit)
    qc.rz(-2 * phi, qubit)
    qc.z(qubit)

    return qc


def insert_error(qc, pos, k):
    """
    Insert independent continuous errors on a set of qubits.

    For every qubit in ``pos``, the two angles are sampled independently
    from Gaussian distributions with zero mean and standard deviation k/3.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to which the errors are appended.
    pos : iterable of int
        Qubit indices receiving an error.
    k : float
        Continuous-noise strength. The Gaussian standard deviation is k/3.

    Returns
    -------
    QuantumCircuit
        The modified noisy circuit.
    """
    for qubit in pos:
        theta = np.random.normal(loc=0.0, scale=k / 3)
        phi = np.random.normal(loc=0.0, scale=k / 3)

        qc = r_error_gate(qc, theta, phi, qubit)

    return qc


# ---------------------------------------------------------------------
# Grover phase oracle
# ---------------------------------------------------------------------

def phase(n, k, target):
    """
    Construct the noisy phase oracle for one marked computational state.

    The oracle applies a phase of -1 to ``target``. The target string is
    reversed because Qiskit uses little-endian qubit ordering.

    The ideal oracle is decomposed into elementary gates. A fresh coherent
    error is then inserted after every single-qubit gate in the decomposed
    circuit.

    Parameters
    ----------
    n : int
        Number of qubits.
    k : float
        Continuous-noise strength.
    target : str
        Marked n-bit computational-basis state.

    Returns
    -------
    QuantumCircuit
        Noisy phase-oracle circuit.
    """
    # Multi-controlled Z gate acting on all n qubits.
    mcz = ZGate().control(n - 1)

    ideal = QuantumCircuit(n)

    # Reverse the bit string to match Qiskit's qubit ordering.
    reversed_target = target[::-1]

    # Map the marked state to |11...1>.
    for qubit, bit in enumerate(reversed_target):
        if bit == "0":
            ideal.x(qubit)

    # Apply a phase of -1 to |11...1>.
    ideal.compose(mcz, range(n), inplace=True)

    # Undo the initial X gates.
    for qubit, bit in enumerate(reversed_target):
        if bit == "0":
            ideal.x(qubit)

    # Decompose the multi-controlled operation into lower-level gates.
    ideal = ideal.decompose(reps=10)

    # Rebuild the circuit while inserting noise after 1-qubit gates.
    noisy = QuantumCircuit(n)

    for instruction, qargs, cargs in ideal.data:
        noisy.append(instruction, qargs, cargs)

        # No fresh noise is inserted after multi-qubit gates.
        if len(qargs) == 1:
            qubit = noisy.qubits.index(qargs[0])
            noisy = insert_error(noisy, [qubit], k)

    return noisy


# ---------------------------------------------------------------------
# Grover diffuser
# ---------------------------------------------------------------------

def diffuser(n, k):
    """
    Construct the noisy Grover diffusion operator.

    The ideal diffuser is

        H^n X^n MCZ X^n H^n.

    After decomposition, a fresh continuous error is inserted after every
    single-qubit gate.

    Parameters
    ----------
    n : int
        Number of qubits.
    k : float
        Continuous-noise strength.

    Returns
    -------
    QuantumCircuit
        Noisy diffusion circuit.
    """
    mcz = ZGate().control(n - 1)

    ideal = QuantumCircuit(n)

    # Map the uniform superposition state to the phase-flip basis.
    ideal.h(range(n))
    ideal.x(range(n))

    # Apply the multi-controlled phase flip.
    ideal.compose(mcz, range(n), inplace=True)

    # Undo the basis transformation.
    ideal.x(range(n))
    ideal.h(range(n))

    # Decompose into elementary operations before adding gate-level noise.
    ideal = ideal.decompose(reps=10)

    noisy = QuantumCircuit(n)

    for instruction, qargs, cargs in ideal.data:
        noisy.append(instruction, qargs, cargs)

        # Insert noise only after single-qubit gates.
        if len(qargs) == 1:
            qubit = noisy.qubits.index(qargs[0])
            noisy = insert_error(noisy, [qubit], k)

    return noisy


# ---------------------------------------------------------------------
# Complete Grover search circuit
# ---------------------------------------------------------------------

def GSA(n, target, k):
    """
    Build a noisy Grover search algorithm circuit.

    The circuit:

    1. prepares the uniform superposition;
    2. applies coherent noise after the initial Hadamard gates;
    3. repeats the phase oracle and diffuser;
    4. measures all qubits.

    Parameters
    ----------
    n : int
        Number of search qubits.
    target : str
        Marked computational-basis state.
    k : float
        Continuous-noise strength.

    Returns
    -------
    QuantumCircuit
        Complete noisy Grover circuit.
    """
    qc = QuantumCircuit(n)

    # Prepare the uniform superposition.
    qc.h(range(n))

    # Add one continuous error after each initial Hadamard gate.
    qc = insert_error(qc, range(n), k)

    qc.barrier()

    # Approximate optimal number of Grover iterations.
    repetitions = int(np.pi / 4 * np.sqrt(2**n) - 0.5)

    for _ in range(repetitions):
        # Apply the noisy marked-state phase oracle.
        qc.compose(phase(n, k, target), range(n), inplace=True)
        qc.barrier()

        # Apply the noisy Grover diffuser.
        qc.compose(diffuser(n, k), range(n), inplace=True)
        qc.barrier()

    qc.measure_all()

    return qc


# ---------------------------------------------------------------------
# Circuit execution
# ---------------------------------------------------------------------

def QuantumSimulator(qc, shots=100):
    """
    Transpile and execute a circuit using the ideal Aer simulator.

    The circuit already contains the sampled continuous-error gates.
    Therefore, no separate Aer noise model is used here.

    Parameters
    ----------
    qc : QuantumCircuit
        Circuit to simulate.
    shots : int, optional
        Number of measurements in this simulator call.

    Returns
    -------
    dict
        Measurement counts indexed by bit string.
    """
    backend = AerSimulator()

    # Let the transpiler select a gate basis supported by Aer.
    tqc = transpile(
        qc,
        backend=backend,
        optimization_level=3,
    )

    # Another possible choice is:
    #
    # tqc = transpile(
    #     qc,
    #     basis_gates=["u", "cx"],
    #     optimization_level=3,
    # )

    job = backend.run(tqc, shots=shots)
    result = job.result()

    return result.get_counts()


def Sum_dict(dict1, dict2):
    """
    Add the values of two count dictionaries.

    Keys missing from one dictionary are treated as having zero counts.

    Parameters
    ----------
    dict1, dict2 : dict
        Measurement-count dictionaries.

    Returns
    -------
    dict
        Merged count dictionary.
    """
    merged_dict = dict1.copy()

    for key, value in dict2.items():
        if key in merged_dict:
            merged_dict[key] += value
        else:
            merged_dict[key] = value

    return merged_dict


def runcircuitnoisy(qc, shots):
    """
    Run the same sampled noisy circuit in repeated batches.

    Each batch contains 100 Aer shots. Therefore, the effective total
    number of measurements is

        100 * shots.

    Important
    ---------
    The circuit is not rebuilt inside this function. Consequently, all
    batches use the same sampled coherent-error realization.

    Parameters
    ----------
    qc : QuantumCircuit
        Previously constructed noisy circuit.
    shots : int
        Number of 100-shot batches.

    Returns
    -------
    dict
        Counts accumulated over all batches.
    """
    counts = {}

    for _ in range(shots):
        batch_counts = QuantumSimulator(qc, shots=100)
        counts = Sum_dict(counts, batch_counts)

    return counts


# ---------------------------------------------------------------------
# Target-state generation
# ---------------------------------------------------------------------

# This alternative function would return all 2**n computational states:
#
# def binary_strings(n: int) -> list[str]:
#     return [format(i, f"0{n}b") for i in range(2**n)]


def binary_strings(n: int) -> list[str]:
    """
    Generate n sampled target strings.

    The function generates one target with exactly k zero bits for each

        k = 1, ..., n.

    Thus, it returns n targets rather than all 2**n possible target states.
    The zero positions are selected randomly.

    For example, for n=3, it returns:

    - one string containing one zero;
    - one string containing two zeros;
    - the string ``000``.

    Parameters
    ----------
    n : int
        Number of bits in each target.

    Returns
    -------
    list[str]
        List containing n sampled target strings.
    """
    targets = []

    for number_of_zeros in range(1, n + 1):
        # Begin with the all-one state.
        bits = ["1"] * n

        # Randomly select positions that will be changed to zero.
        zero_positions = random.sample(
            range(n),
            number_of_zeros,
        )

        for position in zero_positions:
            bits[position] = "0"

        targets.append("".join(bits))

    return targets


# ---------------------------------------------------------------------
# Noise-strength sweep
# ---------------------------------------------------------------------

def RunGSA_AllNoise(a, b, num, shots, n, seed=None):
    """
    Compute the average Grover failure probability over a noise sweep.

    For each noise strength:

    1. generate a noisy Grover circuit for each sampled target;
    2. execute the circuit;
    3. count measurements different from the marked target;
    4. average the failure probability over the sampled targets.

    Parameters
    ----------
    a : float
        Lower exponent of the logarithmic noise range.
    b : float
        Upper exponent of the logarithmic noise range.
    num : int
        Number of logarithmically spaced noise strengths.
    shots : int
        Number of 100-shot batches for each circuit.
    n : int
        Number of qubits.
    seed : int or None, optional
        Seed used for selecting the target bit strings.

        This currently seeds Python's ``random`` module only. It does not
        seed NumPy, which samples the coherent-error angles.

    Returns
    -------
    list[float]
        Average failure probability for each noise strength.
    """
    if seed is not None:
        random.seed(seed)

    # Generate logarithmically spaced continuous-noise strengths.
    p_errors = np.logspace(a, b, num=num)

    failure_probabilities = []

    # Generate n sampled marked states.
    targets = binary_strings(n)

    # Each call made inside runcircuitnoisy uses 100 Aer shots.
    shots_per_call = 100

    for noise_strength in p_errors:
        failures = 0.0
        total_measurements = 0.0

        for target in targets:
            # Constructing GSA here samples a new coherent-error realization
            # for this target and this noise strength.
            noisy_circuit = GSA(
                n=n,
                target=target,
                k=noise_strength,
            )

            counts = runcircuitnoisy(
                noisy_circuit,
                shots,
            )

            # Measurements of the marked state are counted as successful.
            success = counts.get(target, 0)

            total_shots = shots * shots_per_call

            failures += total_shots - success
            total_measurements += total_shots

        # Average failure probability over all selected target states.
        failure_probabilities.append(
            failures / total_measurements
        )

    return failure_probabilities


# ---------------------------------------------------------------------
# Saving results
# ---------------------------------------------------------------------

def save_table(table, filename):
    """
    Append simulation results to a text file.

    Append mode is used because the simulation can be launched several
    times, including through independent parallel jobs. Each execution
    adds a new result row to the same output file.

    Parameters
    ----------
    table : iterable
        Rows to write to the file.
    filename : str
        Output text-file path.
    """
    with open(filename, "a") as file:
        # Note: this header is hard-coded and does not automatically match
        # the ``shots`` value passed to RunGSA_AllNoise.
        file.write("Shots = 100\n")

        for row in table:
            file.write(" ".join(map(str, row)) + "\n")


# ---------------------------------------------------------------------
# Main simulation
# ---------------------------------------------------------------------

# Run a 20-point noise sweep for a six-qubit Grover circuit.
#
# The argument shots=200 means 200 batches of 100 Aer shots, giving
# 20,000 measurements for each target and noise strength.
results = RunGSA_AllNoise(
    a=-3,
    b=-1,
    num=20,
    shots=100,
    n=3,
)

# Append the raw curve to the output file. Append mode allows independent
# parallel executions to accumulate several realizations in one file.
save_table(
    [results],
    "test.txt",
)
