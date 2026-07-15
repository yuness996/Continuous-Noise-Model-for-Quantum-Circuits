"""Grover search under a symmetric stochastic Pauli noise model.

A Pauli channel is applied after the transpiler's one-qubit basis gates. The
script evaluates the average probability of not measuring the marked state.
"""

from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import ZGate
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, pauli_error


def build_grover_circuit(n_qubits: int, target: str) -> QuantumCircuit:
    """Construct an ideal Grover-search circuit for one marked state."""
    if len(target) != n_qubits:
        raise ValueError("The target length must equal the number of qubits.")

    multi_controlled_z = ZGate().control(n_qubits - 1)
    circuit = QuantumCircuit(n_qubits)

    # Prepare the uniform superposition.
    circuit.h(range(n_qubits))
    circuit.barrier()

    # Build the phase oracle. Reverse the string to match Qiskit's
    # little-endian mapping between displayed bitstrings and qubit indices.
    phase_oracle = QuantumCircuit(n_qubits)
    reversed_target = target[::-1]

    for qubit, bit in enumerate(reversed_target):
        if bit == "0":
            phase_oracle.x(qubit)

    phase_oracle.compose(multi_controlled_z, range(n_qubits), inplace=True)

    for qubit, bit in enumerate(reversed_target):
        if bit == "0":
            phase_oracle.x(qubit)

    # Build the standard Grover diffusion operator.
    diffuser = QuantumCircuit(n_qubits)
    diffuser.h(range(n_qubits))
    diffuser.x(range(n_qubits))
    diffuser.compose(multi_controlled_z, range(n_qubits), inplace=True)
    diffuser.x(range(n_qubits))
    diffuser.h(range(n_qubits))

    # Near-optimal integer number of Grover iterations for one marked state.
    repetitions = int(np.pi / 4 * np.sqrt(2**n_qubits) - 0.5)

    for _ in range(repetitions):
        circuit.compose(phase_oracle, range(n_qubits), inplace=True)
        circuit.barrier()
        circuit.compose(diffuser, range(n_qubits), inplace=True)
        circuit.barrier()

    circuit.measure_all()
    return circuit


def all_binary_strings(n_qubits: int) -> List[str]:
    """Return all ``2**n_qubits`` computational-basis bitstrings."""
    return [format(value, f"0{n_qubits}b") for value in range(2**n_qubits)]


def run_pauli_circuit(
    circuit: QuantumCircuit,
    error_probability: float,
    shots: int,
) -> Dict[str, int]:
    """Run a circuit with a symmetric one-qubit Pauli channel.

    The total non-identity probability is ``error_probability`` and is split
    equally between X, Y, and Z. The channel is attached to the one-qubit
    transpiler basis gates ``u1``, ``u2``, and ``u3``.
    """
    if not 0.0 <= error_probability <= 1.0:
        raise ValueError("error_probability must lie between 0 and 1.")

    pauli_channel = pauli_error(
        [
            ("X", error_probability / 3),
            ("Y", error_probability / 3),
            ("Z", error_probability / 3),
            ("I", 1 - error_probability),
        ]
    )

    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(
        pauli_channel,
        ["u1", "u2", "u3"],
    )

    backend = AerSimulator(noise_model=noise_model)
    transpiled = transpile(circuit, backend=backend, optimization_level=3)
    result = backend.run(transpiled, shots=shots).result()
    return result.get_counts()


def run_grover_pauli_sweep(
    exponent_min: float,
    exponent_max: float,
    num_points: int,
    shots: int,
    n_qubits: int,
) -> List[float]:
    """Return the average failure probability versus Pauli error strength.

    For every noise value, all computational-basis targets are tested. The
    failure probability is

        1 - N_target / N_total,

    where ``N_target`` is the number of measurements equal to the marked state.
    """
    error_probabilities = np.logspace(
        exponent_min,
        exponent_max,
        num=num_points,
    )
    failure_probabilities = []
    targets = all_binary_strings(n_qubits)

    for error_probability in error_probabilities:
        failures = 0
        total_measurements = 0

        for target in targets:
            circuit = build_grover_circuit(n_qubits, target)
            counts = run_pauli_circuit(circuit, error_probability, shots)

            successes = counts.get(target, 0)
            failures += shots - successes
            total_measurements += shots

        failure_probabilities.append(failures / total_measurements)

    return failure_probabilities


def plot_results(
    error_probabilities: np.ndarray,
    failure_probabilities: List[float],
) -> None:
    """Plot the Grover failure probability on logarithmic axes."""
    plt.plot(error_probabilities, failure_probabilities, label="Pauli model")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Physical error probability $p$")
    plt.ylabel("Failure probability")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    minimum_exponent = -6
    maximum_exponent = -1
    number_of_points = 20

    results = run_grover_pauli_sweep(
        exponent_min=minimum_exponent,
        exponent_max=maximum_exponent,
        num_points=number_of_points,
        shots=100_000,
        n_qubits=6,
    )

    print(f"GSA_Pauli6=[{', '.join(map(str, results))}]")

    physical_error_probabilities = np.logspace(
        minimum_exponent,
        maximum_exponent,
        num=number_of_points,
    )
    plot_results(physical_error_probabilities, results)